#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Pdbaa;
use Sql;
use Parallel::ForkManager;


use constant ASSEMBLIES => '/home/tmp/BISC/update/pisa/lists/assemblies/';
use constant PDB_STRUCTURES => '/home/tmp/BISC/update/pdb/structures/';
use constant PISA_MODIFIED => '/home/tmp/BISC/update/pisa/structures/modified/';
use constant PDBAA  => '/home/tmp/BISC/update/pdbaa.new';
use constant PISA2PDB => '/home/tmp/BISC/scripts/2010/pisa2pdb.log'; 
use constant PISA2PDB_OUT => '/home/tmp/BISC/scripts/2010/pisa2pdb_DB'; 
use constant PISACFG1 => '/home/tmp/BISC/pisa/mmdb/pisa.bisc1.cfg'; 
use constant PISACFG2 => '/home/tmp/BISC/pisa/mmdb/pisa.bisc2.cfg';


our $sql = new Sql;
our $pdbaa = new Pdbaa;
our $pm = new Parallel::ForkManager(4);

&main();

##########################
# populates pdbPisaChains, everything but pisa protein Chains
##########################

sub main(){
#   my $Lpdbaa = '/home/s0571283/BISC/input/pdbaa.test';
#   print Dumper($HRpisa2pdb);die; 
#my $HRpisa2pdb = &getPisaChains();#print Dumper($HRpisa2pdb);

&populatePdbPisaChains();
}

sub getPisaChains(){
   my ($HRresult);
   my (%log);
   my $ARpisaAssemblies = $pdbaa->parseDir(ASSEMBLIES,'pisa.assemblies$');
#                                                      print Dumper(\@pisa);die; 
   my $HR;
PISA:
   for my $pisa (@{$ARpisaAssemblies}){
      $pm->start and next; 
      my ($session) = (0);

      my @path = split(/\//,$pisa);
      my @last = split(/\./,$path[$#path]); 
      my $pdbID = $last[0];

#get number of protein chains in PDB
      my $LpdbFile = PDB_STRUCTURES . $pdbID . '.pdb';
      if(not -e $LpdbFile){die $LpdbFile}
      my ($HRpdb) = $pdbaa->parseATOMseq($LpdbFile);
      $HRresult->{$pdbID}->{proteinChainsPdb} = keys(%{$HRpdb->{normal}});
      undef($HRpdb);


# If we were able to map chains, a modified file should exist
# Get chains from it.

      my $ARmod = $pdbaa->parseDir(PISA_MODIFIED,"^$pdbID\.pisa\.modified\..*\.pdb");
      if($ARmod->[0]){
         my $Lmod = PISA_MODIFIED . $ARmod->[0];
         my ($HRpdb) = $pdbaa->parseATOMseq($Lmod);
         $HRresult->{$pdbID}->{mappedProteinChainsPisa} = keys(%{$HRpdb->{normal}});
         undef($HRpdb);
         
      }
      else{
         $HRresult->{$pdbID}->{mappedProteinChainsPisa} = 0;
      }


#      next unless($pdbID eq '10gs');
#      if(exists $log{$pdbID}){die $pdbID}
#      else{$log{$pdbID}++}

# information about the assembly, serial number and total chains by PISA      

      my $Lassemblies = ASSEMBLIES . $pdbID . '.pisa.assemblies';
      my $HRassemblies = $pdbaa->parsePisaAssembly($Lassemblies);

#print Dumper($HRassemblies);
#check if session 1 or 2
      if (-e "/home/tmp/BISC/update/pisa/sessions1/pisa_$pdbID"){$session = 1;}
      elsif (-e "/home/tmp/BISC/update/pisa/sessions2/pisa_$pdbID"){$session = 2;}
      else {warn $pdbID;$session = 0;}

# $key indicates if complex is stable. Possibilities: >0,0,-1,-2,-3
# >0 = stable, number indicates the amount of macromolecular chains
# 0 = grey zone, assemblies in this category might be stabl or not
# -1 = unstable assemblies
# -2 = No calculation found in  session
# -3 = Does not complexate in solution
# NULL = no PISA calculation
# if no stable complex exists (1), the next highest number is taken
      if (not exists $HRassemblies->{1}){
         foreach my $chainsPISA (reverse sort keys %{$HRassemblies}){
            $HRresult->{$pdbID}->{totalPISAchains} = $chainsPISA;
            $HRresult->{$pdbID}->{serial} = '0';
            $HRresult->{$pdbID}->{session} = $session;
#print  "          #pdbID | chainsPDB                              | chainsPISA                           | pisaSerial                  | session                        | proteinPisa\n"; 
#            print  "$pdbID,$HRresult->{$pdbID}->{proteinChainsPdb},$HRresult->{$pdbID}->{totalPISAchains},$HRresult->{$pdbID}->{serial},$HRresult->{$pdbID}->{session},$HRresult->{$pdbID}->{mappedProteinChainsPisa}\n ";
#            print  "PDB:$pdbID, PCpdb:$HRresult->{$pdbID}->{proteinChainsPdb}, TotalPisa:$HRresult->{$pdbID}->{totalPISAchains}, Serial:$HRresult->{$pdbID}->{serial}, Session:$HRresult->{$pdbID}->{session}, Mapped:$HRresult->{$pdbID}->{mappedProteinChainsPisa}\n ";
            $HR->{$pdbID} = "$pdbID\t$HRresult->{$pdbID}->{proteinChainsPdb}\t$HRresult->{$pdbID}->{totalPISAchains}\t$HRresult->{$pdbID}->{serial}\t$HRresult->{$pdbID}->{session}\t    $HRresult->{$pdbID}->{mappedProteinChainsPisa}";
            $pm->finish;
            next PISA;
         }

      }
#get the set with the most chains
      else{
         foreach my $size ( reverse sort {$a <=> $b  } keys %{$HRassemblies->{1}->{1} }){
#         print "$size\n"; 
#         print Dumper($HRassemblies->{1}->{1}); 
            my $serial = $HRassemblies->{1}->{1}->{$size}->{Serial};
#            print "$size,$serial,$session,$pdbID\n"; 
            $HRresult->{$pdbID}->{totalPISAchains} = $size;
            $HRresult->{$pdbID}->{serial} = $serial,;
            $HRresult->{$pdbID}->{session} = $session;
            if(not exists $HRresult->{$pdbID}->{mappedProteinChainsPisa}){$HRresult->{$pdbID}->{mappedProteinChainsPisa} = 0}
#print  "pdbID |chainsPDB | chainsPISA | pisaSerial | session | proteinPisa\n"; 
#print  "          #pdbID | chainsPDB                              | chainsPISA                           | pisaSerial                  | session                        | proteinPisa\n"; 
#            print  "PDB:$pdbID, PCpdb:$HRresult->{$pdbID}->{proteinChainsPdb}, TotalPisa:$HRresult->{$pdbID}->{totalPISAchains}, Serial:$HRresult->{$pdbID}->{serial}, Session:$HRresult->{$pdbID}->{session}, Mapped:$HRresult->{$pdbID}->{mappedProteinChainsPisa}\n ";
            $HR->{$pdbID} = "$pdbID\t$HRresult->{$pdbID}->{proteinChainsPdb}\t$HRresult->{$pdbID}->{totalPISAchains}\t$HRresult->{$pdbID}->{serial}\t$HRresult->{$pdbID}->{session}\t    $HRresult->{$pdbID}->{mappedProteinChainsPisa}";
            
            $pm->finish;
#   return($HRresult); 
            next PISA;
         }
      } 
      $pm->finish;
   }
   $pm->wait_all_children;
   return($HR); 
}

sub writePdbPisaChains(\%){
   my ($HR) = @_;
   open(FH,'>',PISA2PDB_OUT) or die "Can not open/access '".PISA2PDB_OUT."'\n$!";
      print FH "#pdbID | chainsPDB | chainsPISA | pisaSerial | session | proteinPisa\n"; 
      foreach my $pdbID ( sort  keys %{$HR}){
         if (not defined $HR->{$pdbID}->{chainsPISA} ){$HR->{$pdbID}->{chainsPISA} = '-1'}
         if (not defined $HR->{$pdbID}->{serial} ){$HR->{$pdbID}->{serial} = '-1'}
         if (not defined $HR->{$pdbID}->{session} ){$HR->{$pdbID}->{session} = '-1'}
         print FH $HR->{$pdbID}; 

      }
   close(FH); 
}

sub populatePdbPisaChains(\%){
   my ($HR) = @_;
   my $dbh;
   $sql->connectBiscTest(\$dbh);
   $dbh->do('truncate pdbPisaChains');
   
   my $sth = $dbh->prepare(q{
         insert into  pdbPisaChains  (pdbID, chainsPDB, chainsPISA, pisaSerial, session, proteinPisa)  values (?,?,?,?,?,?) ;
         }) or die $dbh->errstr;
#pdbID | chainsPDB | chainsPISA | pisaSerial | session | proteinPisa
# 3ihd  1  -3 0  2  0

   open(FH,'<',PISA2PDB_OUT) or die "Can not open/access '".PISA2PDB_OUT."'\n$!";
      while(my $line = <FH>){
         next if($line =~ /^#/);
         chomp($line);
         $line =~ /^(\w{4})\t(\S*)\t(\S*)\t(\S*)\t(\S*)\t(\S*)/;
         $sth->execute($1,$2,$3,$4,$5,$6);
      }
   close(FH);

   $sth->finish();
   $sql->disconnectBisc(\$dbh); 
}
=begin comment
   #my @pisa = </home/tmp/BISC/update/pisa/lists/assemblies/*.pisa.assemblies>;
sub getPdbChains(\%){
   my ($HRpdbaa) = @_;
   my (%size);
#print Dumper($HRpdbaa); 

   foreach my $pdbIDuc (sort keys %{$HRpdbaa->{'experiment'}}){
      next if($HRpdbaa->{'experiment'}->{$pdbIDuc} ne 'XRAY');
         my $chains = scalar(@{$HRpdbaa->{'list'}->{$pdbIDuc}}); 
         my $pdbIDlc = lc($pdbIDuc);
         $size{$pdbIDlc}{chainsPDB} = $chains;
   }
#   print Dumper(\%size); 
#   print "Size: ".keys (%size)."\n";
   return(\%size);
}   
#      if(!$ARmod){
#        if (-e "/home/tmp/BISC/update/pisa/sessions1/pisa_$pdbID" ){
#           my ($fh, $filename) = tempfile();
#           my $args = "pisa $pdbID -list assemblies /home/tmp/BISC/pisa/mmdb/pisa.bisc1.cfg > $filename"; 
#           system($args) == 0 or die $args;
#           open(FH,'<',$filename) or die "Can not open/access '".$filename."'\n$!";
#           while(my $line = <FH>){
#              if($line =~ /quaternary structures are stable in solution/){$stable = '1'}
#              elsif($line =~ /quaternary structure is stable in solution/){$stable = '1'}
#              elsif($line =~ /These structures may or may not be/){$stable = '0'}
#              elsif($line =~ /This structure may or may not be/){$stable = '0'}
#              elsif($line =~ /however they do not/){$stable = '-1'}
#              elsif($line =~ /no calculation results found in session/){$stable = '-2';$result{$stable}{0}{0} = 'no calculation results found in session';next;}
#              elsif($line =~ /complexate in solution/){$stable = '-3';$result{$stable}{0}{0} = 'Does not complexate in solution';next;}
#
#            }
#            close(FH);
#             
#        }
#      }
#      print PISA_MODIFIED. "$pdbID\n"; 
#      print Dumper($ARmod); die; 
