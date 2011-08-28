#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
#use lib '/home/tmp/BISC/libs/';
use Pdbaa;
use Sql;
use Parallel::ForkManager;
use File::Temp qw/ tempfile tempdir /;

our $pdbaa = new Pdbaa;
our $sql = new Sql;
our $pm = new Parallel::ForkManager(4);

use constant DIRECTORIES => '/home/tmp/BISC/update/subcomplexes/';
use constant IFRESIDUES  => '/home/tmp/BISC/update/html/interfaceResidues/';
use constant STRUCTURES  => 'http://bisc.cse.ucsc.edu/structures/'; #This path is has to point to directory where the structures are stored on the webserver
use constant INTERACTINGHOM => '/home/tmp/BISC/update/output/BiscHom/BiscHom.interacting';
use constant INTERACTINGHET => '/home/tmp/BISC/update/output/BiscHet/BiscHet.interacting';
use constant OUTPUTPISA     => '/home/tmp/BISC/update/upload/interfaceResidues/pisa/';
use constant OUTPUTPDB      => '/home/tmp/BISC/update/upload/interfaceResidues/pdb/';
use constant NACCESS        => '/home/tmp/BISC/naccess2.1.1/naccess'; 
use constant PISACFG1       => '/home/tmp/BISC/pisa/mmdb/pisa.bisc1.cfg';
use constant PISACFG2       => '/home/tmp/BISC/pisa/mmdb/pisa.bisc2.cfg';
use constant PDBCFG         => '/home/tmp/BISC/pisa/mmdb/pisa.pdb.cfg';
use constant PISA           => '/home/tmp/BISC/pisa/pisa/pisa';
use constant PISA_SUB_CFG => '/home/tmp/BISC/pisa/mmdb/pisa.subcomplexes.cfg';
use constant PISA_SUB_CFG2 => '/home/tmp/BISC/pisa/mmdb/pisa.subcomplexes.2.cfg';
use constant PISA_SUB_CFG3 => '/home/tmp/BISC/pisa/mmdb/pisa.subcomplexes.3.cfg';

use constant SUB_PDB_HOM   => '/home/tmp/BISC/update/upload/structures/BiscHom/pdb/';
use constant SUB_PISA_HOM  => '/home/tmp/BISC/update/upload/structures/BiscHom/pisa/';
use constant SUB_PDB_HET   => '/home/tmp/BISC/update/upload/structures/BiscHet/pdb/';
use constant SUB_PISA_HET  => '/home/tmp/BISC/update/upload/structures/BiscHet/pisa/';

#use constant SUB_PISA_HOM => '/home/tmp/BISC/update/upload/BiscHom/subcomplexes/pisa/';
#use constant SUB_PDB_HOM => '/home/tmp/BISC/update/upload/BiscHom/subcomplexes/pdb/';
#use constant SUB_PISA_HET => '/home/tmp/BISC/update/upload/BiscHet/subcomplexes/pisa/';
#use constant SUB_PDB_HET => '/home/tmp/BISC/update/upload/BiscHet/subcomplexes/pdb/';


#&openLog();
&main();
#&closeLog();

sub main(){
   my $dbh;
#my $HRlog = &parseLog($Llog);#print Dumper($HRlog); die;
   &hom();
   &het();
}

sub hom(){
   my $HRinteractingHom    = $pdbaa->parseInteractions(INTERACTINGHOM);
   &identifyInterfaceResidues($HRinteractingHom->{'line'},'Hom');
}

sub het(){
   my $HRinteractingHet    = $pdbaa->parseInteractions(INTERACTINGHET);
   &identifyInterfaceResidues($HRinteractingHet->{'line'},'Het');
}


sub identifyInterfaceResidues(\$\%$){ 
   if(scalar @_ != 2){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 2 needed\n";die;}
   my ($HRinteractingFile,$FhomHet) = @_;
   my ($sth,$sthUpdatePdb,$sthUpdatePisa,$subDir,$session);
   my $dbh;
   $sql->connectBisc(\$dbh);
   if($FhomHet eq 'Hom'){
#$dbh->do('truncate biscHom');
      $sth = $dbh->prepare(q{
            insert into biscHom (pdbID, chainID1, chainID2, pdbInterface, pisaInterface, pdbAssembly, pisaAssembly) values (?,?,?,?,?,?,?);
            }) or die $dbh->errstr;
      $sthUpdatePdb = $dbh->prepare(q{
            update biscHom set  pdbDeltaG = ?, pdbNhb = ?, pdbNsb = ?, pdbNds = ?, pdbSymetrie = ?  where pdbID = ? and  chainID1 = ?  and  chainID2 = ?;
            }) or die $dbh->errstr;
      $sthUpdatePisa = $dbh->prepare(q{
            update biscHom set  pisaDeltaG = ?, pisaNhb = ?, pisaNsb = ?, pisaNds = ?, pisaSymetrie = ?  where pdbID = ? and  chainID1 = ?  and  chainID2 = ?;
            }) or die $dbh->errstr;
#
      $subDir = 'BiscHom/';
      $session = 'Hom';
   }
   elsif($FhomHet eq 'Het'){
#$dbh->do('truncate biscHet');
      $sth = $dbh->prepare(q{
            insert into biscHet (pdbID, chainID1, chainID2, pdbInterface, pisaInterface, pdbAssembly, pisaAssembly) values (?,?,?,?,?,?,?);
            }) or die $dbh->errstr;
      $sthUpdatePdb = $dbh->prepare(q{
            update biscHet set  pdbDeltaG = ?, pdbNhb = ?, pdbNsb = ?, pdbNds = ?, pdbSymetrie = ?  where pdbID = ? and  chainID1 = ?  and  chainID2 = ?;
            }) or die $dbh->errstr;
      $sthUpdatePisa = $dbh->prepare(q{
            update biscHet set  pisaDeltaG = ?, pisaNhb = ?, pisaNsb = ?, pisaNds = ?, pisaSymetrie = ?  where pdbID = ? and  chainID1 = ?  and  chainID2 = ?;
            }) or die $dbh->errstr;
      $subDir = 'BiscHet/';
      $session = 'Het';
   }
   else{die $FhomHet}
#   my $HRsessions = &getSessions();
#print Dumper($HRinteractingFile);die; 
   foreach my $pdbID (sort keys %{$HRinteractingFile}){
$pm->start and next; 
      my %var;
      $var{pdbID} = $pdbID;
      $var{dir} = substr($pdbID,1,2) . '/';
      $var{pathPisa} = DIRECTORIES . 'pisa/' . $var{dir};
      $var{pathPdb} = DIRECTORIES . 'pdb/' . $var{dir};

      foreach my $chain1 (sort keys %{$HRinteractingFile->{$pdbID}}){
         $var{chain1} = $chain1;
      
         foreach my $chain2 (sort keys %{$HRinteractingFile->{$pdbID}->{$chain1}}){
            my $HRresult;
            $var{chain2} = $chain2;
            $var{pisa}{Lcomplex}  = $var{pathPisa} . $pdbID . '_' . $chain1 . $chain2 .'.rsa';
            $var{pdb}{Lcomplex}   = $var{pathPdb}  . $pdbID . '_' . $chain1 . $chain2 .'.rsa';
            $var{complexName} = "$pdbID:$chain1:$chain2";
            $HRresult->{pisa}->{interface} = -1;
            $HRresult->{pdb}->{interface} = -1;
            my ($Fpdb,$Fpisa) = (0,0);

            if(-e $var{pdb}{Lcomplex}){
               $Fpdb = 1;
               $var{pdb}{Lsingle1} = $var{pathPdb} . $pdbID . '_' . $chain1  .'.rsa';        
               $var{pdb}{Lsingle2} = $var{pathPdb} . $pdbID . '_' . $chain2  .'.rsa';       
               if(not -e $var{pdb}{Lsingle1}){&runNaccess($var{pdb}{Lsingle1});warn $var{pathPdb} . $pdbID . '_' . $chain1 . '.rsa' . ' does not exist'}
               if(not -e $var{pdb}{Lsingle2}){&runNaccess($var{pdb}{Lsingle2});warn $var{pathPdb} . $pdbID . '_' . $chain2 . '.rsa' . ' does not exist'}
               $HRresult->{pdb} = &findInteractingResidues(\%var,$subDir,'pdb');#print Dumper($HRresult); die;
               
               my $id = $pdbID . "_$chain1" . $chain2. '_Pdb_' . $session;
               my $Lfile; 
               if($session eq 'Hom'){$Lfile = SUB_PDB_HOM . $pdbID . '_' . $chain1 . $chain2 .'.pdb'}
               elsif($session eq 'Het'){$Lfile = SUB_PDB_HET . $pdbID . '_' . $chain1 . $chain2 .'.pdb'}
               my $tmp = File::Temp->new( 
                     TEMPLATE => 'tempXXXXX',
                     SUFFIX   => '.tmp');
               my $tempPisa = $tmp->filename;
               &runPisa($id,$Lfile,$tempPisa);
               &parsePisaInterface($tempPisa,$HRresult->{pdb}->{interface},'pdb',\%var,$HRresult->{pdb});
               my $args = "rm -f $tempPisa";
               system ($args) == 0 or warn "$args";
#               print Dumper($HRresult); 
            }

            if(-e $var{pisa}{Lcomplex}){
               $Fpisa = 1;
               $var{pisa}{Lsingle1} = $var{pathPisa} . $pdbID . '_' . $chain1  .'.rsa';        
               $var{pisa}{Lsingle2} = $var{pathPisa} . $pdbID . '_' . $chain2  .'.rsa';       
               if(not -e $var{pisa}{Lsingle1}){&runNaccess($var{pisa}{Lsingle1});warn $var{pathPisa} . $pdbID . '_' . $chain1 . '.rsa' . ' does not exist'}
               if(not -e $var{pisa}{Lsingle2}){&runNaccess($var{pisa}{Lsingle2});warn $var{pathPisa} . $pdbID . '_' . $chain2 . '.rsa' . ' does not exist'}
               $HRresult->{pisa} = &findInteractingResidues(\%var,$subDir,'pisa');#print Dumper($HRresult); die;

               my $id = $pdbID . "_$chain1" . $chain2. '_Pisa_' . $session;
               my $Lfile; 
               if($session eq 'Hom'){$Lfile    = SUB_PISA_HOM . $pdbID . '_' . $chain1 . $chain2 .'.pdb'}
               elsif($session eq 'Het'){$Lfile = SUB_PISA_HET . $pdbID . '_' . $chain1 . $chain2 .'.pdb'}
               my $tmp = File::Temp->new( 
                     TEMPLATE => 'tempXXXXX',
                     SUFFIX => '.tmp');
               my $tempPisa = $tmp->filename;
               &runPisa($id,$Lfile,$tempPisa);
               &parsePisaInterface($tempPisa,$HRresult->{pisa}->{interface},'pisa',\%var,$HRresult->{pisa});
               my $args = "rm -f $tempPisa";
               system ($args) == 0 or warn "$args";
#               print Dumper($HRresult); 
            }
            if( ( $HRresult->{pdb}->{interface} < 200 )  && ( $HRresult->{pisa}->{interface} > 200 )  ){
               &printOutput($HRresult->{pisa},OUTPUTPISA,\%var,'pisa');
            }
            elsif( ( $HRresult->{pdb}->{interface} < 200 )  && ( $HRresult->{pisa}->{interface} > 200 )  ){
               &printOutput($HRresult->{pisa},OUTPUTPISA,\%var,'pisa');
            }
            elsif( ( $HRresult->{pdb}->{interface} > 200 )  && ( $HRresult->{pisa}->{interface} < 200 )  ){
               &printOutput($HRresult->{pdb},OUTPUTPDB,\%var,'pdb');
            }
            elsif( ( $HRresult->{pdb}->{interface} == -1 )  ||  ( $HRresult->{pdb}->{interface} == 0 )  ){
               &printOutput($HRresult->{pisa},OUTPUTPISA,\%var,'pisa');
            }
            elsif( ( $HRresult->{pisa}->{interface} == -1 )  ||  ( $HRresult->{pisa}->{interface} == 0 )  ){
               &printOutput($HRresult->{pdb},OUTPUTPDB,\%var,'pdb');
            }
            else{
               &printOutput($HRresult->{pisa},OUTPUTPISA,\%var,'pisa');
            }
=begin comment
#            print "$pdbID,$chain1,$chain2,\nPDB: $HRresult->{pdb}->{interface},\nPISA:$HRresult->{pisa}->{interface}\n"; 
               $sth->execute($pdbID,$chain1,$chain2,$HRresult->{pdb}->{interface},$HRresult->{pisa}->{interface},$Fpdb,$Fpisa) or die  $dbh->errstr;
            
            if(-e $var{pdb}{Lcomplex}){
               $sthUpdatePdb->execute($HRresult->{pdb}->{deltaG},$HRresult->{pdb}->{Nhb},$HRresult->{pdb}->{Nsb},$HRresult->{pdb}->{Nds},$HRresult->{pdb}->{symOp},$pdbID,$chain1,$chain2);
            }
            if(-e $var{pisa}{Lcomplex}){
               $sthUpdatePisa->execute($HRresult->{pisa}->{deltaG},$HRresult->{pisa}->{Nhb},$HRresult->{pisa}->{Nsb},$HRresult->{pisa}->{Nds},$HRresult->{pisa}->{symOp},$pdbID,$chain1,$chain2);
            }
=cut            
         }   
      }
$pm->finish;
   }
$pm->wait_all_children;
   $sql->disconnectBisc(\$dbh); 
}



#print Dumper($HRresult);sleep(1); 
# opens complex (1FNT_AB.rsa) and single files (1FNT_A, 1FNT_B)
# ASA of each residue in complex file is being compared with each in the single files
# If it differs it is considered part of the interface
sub findInteractingResidues($){
   if(scalar @_ != 3){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 3 needed\n";die;}

   my ($HRvar,$subDir,$FpdbPisa) = @_;
   my (%result,$HRrsa,$args,$path,$fileOnServer);
   my $interface = -1;

   if ($FpdbPisa eq 'pdb'){
      $fileOnServer =  STRUCTURES . $subDir . 'pdb/'  . $HRvar->{pdbID} . '_' . $HRvar->{chain1} . $HRvar->{chain2} .'.pdb';
      $HRrsa->{complex} = $pdbaa->parseNaccessRsa($HRvar->{pdb}->{Lcomplex}); 
      $HRrsa->{single1} = $pdbaa->parseNaccessRsa($HRvar->{pdb}->{Lsingle1});
      $HRrsa->{single2} = $pdbaa->parseNaccessRsa($HRvar->{pdb}->{Lsingle2});
      $result{interface} = &calculateInterface($HRrsa,$HRvar,'pdb');
      if ($result{interface} == 0){
         $result{jmol1}{$HRvar->{complexName}}  =  '"background white; load '. $fileOnServer .'; spacefill off; wireframe off; ';
         $result{'list'}{$HRvar->{complexName}} = -1;
         return(\%result);
      }
   }
   elsif ($FpdbPisa eq 'pisa'){
      $fileOnServer =  STRUCTURES . $subDir . 'pisa/'  . $HRvar->{pdbID} . '_' . $HRvar->{chain1} . $HRvar->{chain2} .'.pdb';
      $HRrsa->{complex} = $pdbaa->parseNaccessRsa($HRvar->{pisa}->{Lcomplex}); 
      $HRrsa->{single1} = $pdbaa->parseNaccessRsa($HRvar->{pisa}->{Lsingle1});
      $HRrsa->{single2} = $pdbaa->parseNaccessRsa($HRvar->{pisa}->{Lsingle2});
      $result{interface} = &calculateInterface($HRrsa,$HRvar,'pisa');
      if ($result{interface} == 0){
         $result{jmol1}{$HRvar->{complexName}}  =  '"background white; load '. $fileOnServer .'; spacefill off; wireframe off; ';
         $result{'list'}{$HRvar->{complexName}} = -1;
         return(\%result);
      }
   }
   else{die $FpdbPisa}
      

#if($complex eq '1b25:A:D'){print "\nC:  $fileC\nS1: $fileS1\nS2: $fileS2\n$path\n$structures\n"; }
=begin comment
   jmolApplet([frameWidth* .40,frameHeight* .80],"background white; load ../../structures/BiscHom/1PMA_IJ.pdb; spacefill off; wireframe off;  select *:I; color [xFCDFFF]; cartoon on; select *:J; color [xECE5B6];  cartoon on;select LYS114:I, LYS41:I; color [x15317E]; wireframe 50 on; spacefill 80;select GLN119:J,ARG93:J;color [xFF0000]; wireframe 50 on; spacefill 80;  select protein; cartoon on; color cartoons chain;");

 jmolApplet([frameWidth* .40,frameHeight* .80],"background white; load ../../structures/BiscHom/1PMA_IJ.pdb; spacefill off; wireframe off;  select *:I; isosurface i1 saSurface; color isosurface  [xFCDFFF]; color isosurface translucent; select *:J; isosurface i2 saSurface; color isosurface [xECE5B6]; color isosurface translucent;select LYS114:I, LYS41:I;color [x15317E]; spacefill 500;select GLN119:J, ARG93:J;color [xFF0000]; spacefill 500;");
=cut 

   $result{jmol1}{$HRvar->{complexName}}  = 'jmolApplet([frameWidth* .40,frameHeight* .80],"background white; load '. $fileOnServer .'; spacefill off; wireframe off; ';
   $result{jmol2}{$HRvar->{complexName}}  = 'jmolApplet([frameWidth* .40,frameHeight* .80],"background white; load '. $fileOnServer .'; spacefill off; wireframe off; ';
   $result{jmol1}{$HRvar->{complexName}} .= " select *:$HRvar->{chain1}; color [xFCDFFF]; cartoon on; select *:$HRvar->{chain2}; color [xECE5B6];  cartoon on;";
   $result{jmol2}{$HRvar->{complexName}} .= " select *:$HRvar->{chain1}; isosurface i1 saSurface; color isosurface [xFCDFFF]; color isosurface translucent; select *:$HRvar->{chain2}; isosurface i2 saSurface; color isosurface [xECE5B6]; color isosurface translucent;";
   $result{outjmol1}{$HRvar->{complexName}}  = "background white;\nload $fileOnServer;\nspacefill off;\nwireframe off;\n";
   $result{outjmol2}{$HRvar->{complexName}}  = "background white;\nload $fileOnServer;\nspacefill off;\nwireframe off;\n";
   $result{outjmol1}{$HRvar->{complexName}} .= "select *:$HRvar->{chain1};\ncolor [xFCDFFF];\ncartoon on;\nselect *:$HRvar->{chain2};\ncolor [xECE5B6];\ncartoon on;\n";
   $result{outjmol2}{$HRvar->{complexName}} .= "select *:$HRvar->{chain1};\nisosurface i1 saSurface;\ncolor isosurface [xFCDFFF];\ncolor isosurface translucent;\nselect *:$HRvar->{chain2};\nisosurface i2 saSurface;\ncolor isosurface [xECE5B6];\ncolor isosurface translucent;\n";

  
# find interface  residues of chain1, add select to Jmol command
   &compareSingleComplex($HRvar->{complexName},$HRrsa->{complex},$HRrsa->{single1},$HRvar->{chain1},\%result,$HRvar->{Lcomplex});
   $result{jmol1}{$HRvar->{complexName}} .= " color [x15317E]; wireframe 50 on; spacefill 80; "; 
   $result{jmol2}{$HRvar->{complexName}} .= " color [x15317E]; spacefill 500; "; 
   $result{outjmol1}{$HRvar->{complexName}} .= "color [x15317E];\nwireframe 50 on;\nspacefill 80;\n"; 
   $result{outjmol2}{$HRvar->{complexName}} .= "color [x15317E];\nspacefill 500;\n"; 
  
# find interface  residues of chain2 
   &compareSingleComplex($HRvar->{complexName},$HRrsa->{complex},$HRrsa->{single2},$HRvar->{chain2},\%result, $HRvar->{Lcomplex});
   $result{jmol1}{$HRvar->{complexName}} .= " color [xFF0000]; wireframe 50 on; spacefill 80; "; 
   $result{jmol2}{$HRvar->{complexName}} .= ' color [xFF0000]; spacefill 500;")'; 
   $result{outjmol1}{$HRvar->{complexName}} .= "color [xFF0000];\nwireframe 50 on;\nspacefill 80;\n"; 
   $result{outjmol2}{$HRvar->{complexName}} .= "color [xFF0000];\nspacefill 500;\n"; 

   $result{jmol1}{$HRvar->{complexName}} .= ' wireframe on; spacefill 80; select protein; cartoon on; color cartoons chain;")';
   $result{outjmol1}{$HRvar->{complexName}} .= "wireframe on;\nspacefill 80;\nselect protein;\ncartoon on;\ncolor cartoons chain;\n";
   
   return(\%result);
}

sub compareSingleComplex(\%$$){  
   my ($complex,$HRrsaComplex,$HRrsaSingle,$chain,$HRresult,$fileC)  = @_;
#      if($complex eq '1b25:A:D'){print "In\n"; } 
   $HRresult->{jmol1}->{$complex} .= " select ";
   $HRresult->{jmol2}->{$complex} .= " select ";
   $HRresult->{outjmol1}->{$complex} .= "select ";
   $HRresult->{outjmol2}->{$complex} .= "select ";
   foreach my $resNo (sort{$a <=> $b} keys %{$HRrsaComplex->{'RES'}->{$chain}}){
#      if($complex eq '1b25:A:D'){print "$resNo\n"; } 
      foreach my $resName (keys %{$HRrsaComplex->{'RES'}->{$chain}->{$resNo}}){
         my $asaC  = $HRrsaComplex->{'RES'}->{$chain}->{$resNo}->{$resName};

#         if($complex eq '1b25:A:D'){print "$resName\t$asaC\n"; } 

         if (defined $HRrsaSingle->{'RES'}->{$chain}->{$resNo}->{$resName}){
            my $asaS = $HRrsaSingle->{'RES'}->{$chain}->{$resNo}->{$resName};
            if($asaS != $asaC){
               my $size = $asaS - $asaC;
               $HRresult->{jmol1}->{$complex} .= $resName . $resNo .":$chain, ";
               $HRresult->{jmol2}->{$complex} .= $resName . $resNo .":$chain, ";
               $HRresult->{outjmol1}->{$complex} .= $resName . $resNo .":$chain, ";
               $HRresult->{outjmol2}->{$complex} .= $resName . $resNo .":$chain, ";
               $HRresult->{'list'}->{$complex} .= "$chain\t$resNo\t$resName\t$size\t$asaS\t$asaC\n";
               $HRresult->{'outlist'}->{$complex} .= "$chain\t$resNo\t$resName\n";

            }
         }
         else{warn "$fileC:$chain:$resNo:$resName does not exist in either single chain file"}
      }
   } 
   $HRresult->{jmol1}->{$complex} =~ s/, $/;/;
   $HRresult->{jmol2}->{$complex} =~ s/, $/;/;
   $HRresult->{outjmol1}->{$complex} =~ s/, $/;/;
   $HRresult->{outjmol2}->{$complex} =~ s/, $/;/;
   $HRresult->{outjmol1}->{$complex} .= "\n";
   $HRresult->{outjmol2}->{$complex} .= "\n";
#die if ($complex eq '1b25:A:D');
}

sub printOutput($){
   if(scalar @_ != 4){die "Argument list wrong: 4"}

   my ($HRresults,$dir,$HRvar,$pdbPisa) =  @_;

   foreach my $id (sort keys %{$HRresults->{jmol1}}){
      if(not defined $HRresults->{'list'}->{$id} ){
         print Dumper($HRvar);
         print "$pdbPisa\n"; die; 
      }

      my @ids = split(/:/,$id);
      my $filename1 = $dir . "$ids[0]_$ids[1]" . $ids[2] . '.select';
      my $filename2 = $dir . "$ids[0]_$ids[1]" . $ids[2] . '.select2';
      my $filename3 = $dir . "$ids[0]_$ids[1]" . $ids[2] . '.wireframe.jmol';
      my $filename4 = $dir . "$ids[0]_$ids[1]" . $ids[2] . '.isosurface.jmol';
      my $filename5 = $dir . "$ids[0]_$ids[1]" . $ids[2] . '.ifResidues';
      my $filename6 = $dir . "$ids[0]_$ids[1]" . $ids[2] . '.interfaceResidues';
      open(FH,'>', $filename1) or die "Can not open/access '" . $filename1 . "'\n$!\n";
         print FH "$HRresults->{jmol1}->{$id}"; 
      close(FH);
      open(FH,'>', $filename2) or die "Can not open/access '" . $filename2 . "'\n$!\n";
         print FH "$HRresults->{jmol2}->{$id}"; 
      close(FH);
      open(FH,'>', $filename3) or die "Can not open/access '" . $filename3 . "'\n$!\n";
         print FH "#Copy and paste the commands in this file in the command window of your local Jmol application.\n"; 
         print FH "$HRresults->{outjmol1}->{$id}"; 
      close(FH);
      open(FH,'>', $filename4) or die "Can not open/access '" . $filename4 . "'\n$!\n";
         print FH "#Copy and paste the commands in this file in the command window of your local Jmol application.\n"; 
         print FH "$HRresults->{outjmol2}->{$id}"; 
      close(FH);
#no interface residues       
      if($HRresults->{'list'}->{$id} eq '-1'){
         warn "$id\nNo interface residues\n";
         return;
      }

      open(FH,'>', $filename5) or die "Can not open/access '" . $filename5 . "'\n$!\n";
         print FH "#Chain | ResidueNumber |  ResidueName | IFsize |  ASAsingle |  ASAcomplex\n"; 
         print FH "$HRresults->{'list'}->{$id}"; 
      close(FH);
      $filename1 =~ s/\.ifResidues/\.interfaceResidues/;

      open(FH,'>', $filename6) or die "Can not open/access '" . $filename6 . "'\n$!\n";
         print FH "#Chain | ResidueNumber |  ResidueName\n"; 
         print FH "$HRresults->{'outlist'}->{$id}"; 
      close(FH);
   } 
}
sub openLog(){
   my $stamp = $pdbaa->stamp();
   print "$stamp\n";
   print STDERR "$stamp\n";

   open(LOG,'>>',$0 . '.log') or die "Can not open/access '".$0.".log'\n$!";
   select((select(LOG), $|=1)[0]); 
   print LOG "$stamp\n";
   
}
sub closeLog(){
   close(LOG);
}
sub runNaccess(){
   my ($path) = @_;
   $path =~ /(.*)_\w{1,2}\.rsa/;
   my $LpdbFile = $1;
   $LpdbFile .= '.pdb';  
   $path =~ /(.*)\w{4}_\w{1,2}\.rsa/;
   my $dir = $1;
   chdir($dir);

   my ($HRpdbFile) = $pdbaa->parsePDBatom($LpdbFile);#print Dumper($HRpdbFile); 

   $path =~ /_(\w{1,2})/;
   if(length($1) == 1){
      my $chain = $1;
      my $Lsingle =  $path;
      $Lsingle =~ s/rsa$/pdb/;

      open(FH,'>',$Lsingle) or die "Can not open/access '$Lsingle'\n$!";
         print FH $HRpdbFile->{'0'}->{$chain};
         print FH "END";
      close(FH);
      my $args = NACCESS . " $Lsingle";
      print "NACCESS: $args\n"; 
      system($args) == 0 or die "Could not $args : $?";
   }
   elsif(length($1) == 2){
      my @chains = split(//,$1);
      my $Lcomplex = $path;
      $Lcomplex =~ s/rsa$/pdb/;
      open(FH,'>',$Lcomplex) or die "Can not open/access '$Lcomplex'\n$!";
         print FH $HRpdbFile->{'0'}->{$chains[0]};
         print FH "TER\n"; 
         print FH $HRpdbFile->{'0'}->{$chains[1]};
         print FH "END";
      close(FH);

      my $args = NACCESS . " $Lcomplex";
      print "NACCESS: $args\n"; 
      system($args) == 0 or die "Could not $args : $?";
   }



}
sub calculateInterface(){
   my ($HRrsa,$HRvar,$FpdbPisa) = @_;
   my ($asaComplex1,$asaComplex2,$asaSingle1,$asaSingle2);
   
   $asaComplex1 = $HRrsa->{complex}->{CHAIN}->{$HRvar->{chain1}}; 
   if( !$asaComplex1 ){ 
      if($FpdbPisa eq 'pdb'){
#         my $args = "ls -l $HRvar->{pathPdb}". $HRvar->{pdbID} . '*';
#         system($args) == 0 or die $args;
         &runNaccess($HRvar->{pdb}->{Lcomplex});
#         system($args) == 0 or die $args;
      }
      elsif($FpdbPisa eq 'pisa'){
#         my $args = "ls -l $HRvar->{pathPisa}". $HRvar->{pdbID} . '*';
#         system($args) == 0 or die $args;
         &runNaccess($HRvar->{pisa}->{Lcomplex});
#         system($args) == 0 or die $args;
      }
      else{die $FpdbPisa}
      $asaComplex1 = $HRrsa->{complex}->{CHAIN}->{$HRvar->{chain1}}; 
   }
   
   $asaComplex2 = $HRrsa->{complex}->{CHAIN}->{$HRvar->{chain2}}; 
   if( !$asaComplex2 ){ 
      if($FpdbPisa eq 'pdb'){
#         my $args = "ls -l $HRvar->{pathPdb}". $HRvar->{pdbID} . '*';
#         system($args) == 0 or die $args;
         &runNaccess($HRvar->{pdb}->{Lcomplex});
#         system($args) == 0 or die $args;
      }
      elsif($FpdbPisa eq 'pisa'){
#         my $args = "ls -l $HRvar->{pathPisa}". $HRvar->{pdbID} . '*';
#         system($args) == 0 or die $args;
         &runNaccess($HRvar->{pisa}->{Lcomplex});
#         system($args) == 0 or die $args;
      }
      else{die $FpdbPisa}
      $asaComplex2 = $HRrsa->{complex}->{CHAIN}->{$HRvar->{chain2}}; 
   }
   
   $asaSingle1  = $HRrsa->{single1}->{CHAIN}->{$HRvar->{chain1}}; 
   if( !$asaSingle1 ){ 
      if($FpdbPisa eq 'pdb'){
#         my $args = "ls -l $HRvar->{pathPdb}". $HRvar->{pdbID} . '*';
#         system($args) == 0 or die $args;
         &runNaccess($HRvar->{pdb}->{Lsingle1});
#         system($args) == 0 or die $args;
      }
      elsif($FpdbPisa eq 'pisa'){
#         my $args = "ls -l $HRvar->{pathPisa}". $HRvar->{pdbID} . '*';
#         system($args) == 0 or die $args;
         &runNaccess($HRvar->{pisa}->{Lsingle1});
#         system($args) == 0 or die $args;
      }
      else{die $FpdbPisa}
      $asaSingle1  = $HRrsa->{single1}->{CHAIN}->{$HRvar->{chain1}}; 
   }

   $asaSingle2  = $HRrsa->{single2}->{CHAIN}->{$HRvar->{chain2}}; 
   if( !$asaSingle2 ){ 
      if($FpdbPisa eq 'pdb'){
#         my $args = "ls -l $HRvar->{pathPdb}". $HRvar->{pdbID} . '*';
#         system($args) == 0 or die $args;
         &runNaccess($HRvar->{pdb}->{Lsingle2});
#        system($args) == 0 or die $args;
      }
      elsif($FpdbPisa eq 'pisa'){
#         my $args = "ls -l $HRvar->{pathPisa}". $HRvar->{pdbID} . '*';
#         system($args) == 0 or die $args;
         &runNaccess($HRvar->{pisa}->{Lsingle2});
#         system($args) == 0 or die $args;
      }
      else{die $FpdbPisa}
      $asaSingle2  = $HRrsa->{single2}->{CHAIN}->{$HRvar->{chain2}}; 
   }

   my $diff1 = $asaSingle1 - $asaComplex1;
   my $diff2 = $asaSingle2 - $asaComplex2;
   my $interface = ($diff1 + $diff2) / 2 ; 
   $interface = sprintf("%d",$interface);
   return($interface);
}
sub parseLog($){
   my $Lfile = shift @_;
   my %result;
   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
   while(my $line = <FH>){
      if($line =~ /^Printed:\t<(.*)>/){
         my @ids = split(/\t/,$1);
         $result{$ids[0]}{$ids[1]}{$ids[2]}++;
      }
   }
   close(FH);

   return(\%result);
}


sub parsePisaInterface(){
   my ($Lfile,$interface,$pdbPisa,$HRvar,$HRresult) = @_;
   my ($difference);

   my ($pdbIDFile);
      my $tmp;
   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
   while(my $line = <FH>){
      if($line =~ /no calculation results/){
         $HRresult->{'areaPisa'} = 'NULL';
         $HRresult->{'areaPdb'}  = 'NULL';
         $HRresult->{'deltaG'}   = 'NULL';
         $HRresult->{'Nhb'}      = 'NULL';
         $HRresult->{'Nsb'}      = 'NULL';
         $HRresult->{'Nds'}      = 'NULL';
         $HRresult->{'symOp'}    = 'NULL';
         return();
      }
      if($line =~ /Session\s(\w{4})/){
         $pdbIDFile = $1; 
         if($HRvar->{pdbID} ne $pdbIDFile){
            die "$HRvar->{pdbID}\t$pdbIDFile"
         }
      }
      next unless ($line =~ /^\s*(\d+)\s+(\d+)\s*\|\s*(\w)[^|]+\|\s*(\w)\s*(\S+)[^|]+\|\s*(\S+)\s*(\S+)\s*(\S+)\s*(\d*)\s*(\d*)\s*(\d*)/ );
      my @chains;
      my $serial = $1;
      push(@chains,$3);
      push(@chains,$4);
      my $sym   = $5; 
      my $area  = $6;
      my $dG    = $7;
      my $nhb   = $8;
      my $nsb   = $9;
      my $nds    = $10;
# next if($area <200);

      @chains = sort(@chains);
      next if( ($HRvar->{chain1} ne $chains[0]) || ($HRvar->{chain2} ne $chains[1]) );
         next unless( $sym eq 'X,Y,Z');

         $area =~ s/\.\d+$//;
         $HRresult->{'areaPisa'} = $area;
         $HRresult->{'areaPdb'} = $interface;
         $HRresult->{'deltaG'} = $dG;
         $HRresult->{'Nhb'} = $nhb;
         $HRresult->{'Nsb'} = $nsb;
         $HRresult->{'Nds'} = $nds;
         $HRresult->{'symOp'} = $sym;

   }
   close(FH);
}

sub runPisa(){
   if(scalar @_ != 3){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 3 needed\n";die;}
   my ($session, $Lfile,$tempPisa) = @_;
   my $args;
#      $args = PISA . " $session -analyse $Lfile " . PISA_SUB_CFG3;
#      print "$args\n"; 
   if(   (not -d '/home/tmp/BISC/update/sessions/subcomplexes/pisa_' . $session) &&  
         (not -d '/home/tmp/BISC/update/sessions/subcomplexes2/pisa_' . $session) && 
         (not -d '/home/tmp/BISC/update/sessions/subcomplexes3/pisa_' . $session)  
         ){
      $args = PISA . " $session -analyse $Lfile " . PISA_SUB_CFG2;
      system($args) == 0 or warn "$args";
   }
#   else{print "$session exists\n"; }
#   my $tmp_fh = new File::Temp( UNLINK => 1 );
#   print "$tmp_fh\n"; 
#     $args = PISA . " $session -list interfaces " . $cfg . " >$tempPisa";
#      system($args) == 0 or warn "$args";
   if( -d '/home/tmp/BISC/update/sessions/subcomplexes/pisa_' . $session  ){
      $args = PISA . " $session -list interfaces " . PISA_SUB_CFG . " >$tempPisa";
      system($args) == 0 or warn "$args";
   }
   elsif( -d '/home/tmp/BISC/update/sessions/subcomplexes2/pisa_' . $session  ){
      $args = PISA . " $session -list interfaces " . PISA_SUB_CFG2 . " >$tempPisa";
      system($args) == 0 or warn "$args";
   }
   elsif( -d '/home/tmp/BISC/update/sessions/subcomplexes3/pisa_' . $session  ){
      $args = PISA . " $session -list interfaces " . PISA_SUB_CFG3 . " >$tempPisa";
      system($args) == 0 or warn "$args";
   }
   else{die PISA . " $session -analyse $Lfile " . PISA_SUB_CFG3 .'did not work'}
   return(1);
}


=cut
'A' => {
   '559' => {
      'ILE' => '1.29'
   },
      '127' => {
         'PHE' => '22.13'
      },
      '32' => {
         'ARG' => '34.20'
      },
      '443' => {
         'VAL' => '6.27'
      },
      '206' => {
         'ALA' => '3.51'
      },
      '118' => {
         'ASN' => '104.51'
      },
      '71' => {
         'ASP' => '119.47'
      },
      '358' => {
         'SER' => '102.01'
      },
      '331' => {
         'CYS' => '0.00'
      },
      '560' => {
         'GLN' => '18.95'
      },

