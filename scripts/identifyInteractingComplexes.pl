#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Pod::Usage;
use Pdbaa;
use Sql;
use Parallel::ForkManager;
use XML::LibXML;
use Carp;


pod2usage if ($ARGV[0] && ($ARGV[0] eq '-h' or $ARGV[0] eq '--help'));

use constant PDBAA           => '/home/tmp/BISC/update/pdbaa.new';
use constant WATER    => '/home/tmp/BISC/EMBOSS-6.0.1/emboss/water';
use constant STARTDIR_PISA  => '/home/tmp/BISC/update/subcomplexes/pisa/';
use constant REFCHAINS => '/home/tmp/BISC/scripts/pdb2pisaChains.log';

use constant PISA_OUTPUT_HOM => '/home/tmp/BISC/update/output/BiscHom/PisaPreCullBiscHom';
use constant PISA_OUTPUT_HET => '/home/tmp/BISC/update/output/BiscHet/PisaPreCullBiscHet';

our $pdbaa = new Pdbaa;
our $pm = new Parallel::ForkManager(4);

my %superseeded = ('2x1y','2xap', '2ehi','3ah1', '2ixw','2xe2', '2ixx','2xe1', '2qyx','3nek');

my $xml_parser  = XML::LibXML->new();
my $cfg         = $xml_parser->parse_file('../config.xml');
#&pisa(\%superseeded);
&pdb(\%superseeded, $cfg);

sub pisa(){
   my ($HRsuperseeded) = @_;
   my $HRrefChains = &parsePisa2PdbMapping(REFCHAINS);
#   print Dumper($HRrefChains); die;
   my $HRpdbaa = $pdbaa->parseDunbrackFasta(PDBAA);
   my $AR = $pdbaa->parseDir(STARTDIR_PISA,'\w\w');

   my $HRoldRun = &parsePisaOutput();
   open(FH,'>',PISA_OUTPUT_HOM) or die "Can not open/access '".PISA_OUTPUT_HOM."'\n$!";
      print FH "#PDB ID | \%seq id | alignment length | indentities  |  ASA interface  | ASA single chain 1  |  ASA complex chain 1 | ASA single chain 2  |  ASA complex chain 2  \n";
   close(FH);
   open(FH,'>',PISA_OUTPUT_HET) or die "Can not open/access '".PISA_OUTPUT_HET."'\n$!";
      print FH "#PDB ID | \%seq id | alignment length | indentities  |  ASA interface  | ASA single chain 1  |  ASA complex chain 1 | ASA single chain 2  |  ASA complex chain 2  \n";
   close(FH);

   for my $sdir (@{$AR}){
      my $path = STARTDIR_PISA . "$sdir/";
      my $ARfiles = $pdbaa->parseDir($path,'\w{4}_\w\w\.rsa');
      for my $rsaFile (@{$ARfiles}){
#      $pm->start and next;
#         if(exists $HRoldRun->{$rsaFile}){
#            $pm->finish;
#            next;

#         }
         $rsaFile =~ /^(\w{4})_(\w)(\w)/;
         my $pdbID = $1;
         my $chain1 = $2;
         my $chain2 = $3;
         next if(exists $HRsuperseeded->{$pdbID});

         my $Lcomplex = $path . $rsaFile;
         my $Lsingle1 = $path . $pdbID . "_$chain1.rsa";
         my $Lsingle2 = $path . $pdbID . "_$chain2.rsa";

         if(not -e $Lcomplex){warn "No\t$Lcomplex";next}
         if(not -e $Lsingle1){warn "No\t$Lsingle1";next}
         if(not -e $Lsingle2){warn "No\t$Lsingle2";next}

         my $asaChain1Complex = &parseNaccessRsa($Lcomplex,$chain1);
         if(! $asaChain1Complex){print STDERR "NACCESS: No score single $pdbID:$chain1\t$Lcomplex\n"; next;}

         my $asaChain2Complex = &parseNaccessRsa($Lcomplex,$chain2);
         if(! $asaChain2Complex){print STDERR "NACCESS: No score single $pdbID:$chain2\t$Lcomplex\n"; next;}

         my $asaChain1Single = &parseNaccessRsa($Lsingle1,$chain1);
         if(! $asaChain1Single){print STDERR "NACCESS: No score single $pdbID:$chain1\t$Lsingle2\n"; next;}

         my $asaChain2Single = &parseNaccessRsa($Lsingle2,$chain2);
         if(! $asaChain2Single){print STDERR "NACCESS: No score single $pdbID:$chain2\t$Lsingle2\n"; next;}

         my $diff1 = $asaChain1Single - $asaChain1Complex;
         my $diff2 = $asaChain2Single - $asaChain2Complex;

         my $interface = ($diff1 + $diff2) / 2 ;
            $interface = sprintf("%d",$interface);
         if($interface > 0){
            my $pdbIDuc = uc($pdbID);
            my ($Fhom,$identities, $alignmentLength, $id) = &classifyInteraction($HRpdbaa,$pdbIDuc,$chain1,$chain2,$HRrefChains);
            if($Fhom == -1){next}

            if($Fhom == 1){
               open(FH,'>>',PISA_OUTPUT_HOM) or die "Can not open/access '".PISA_OUTPUT_HOM."'\n$!";
               select((select(FH), $|=1)[0]);
               print FH "$pdbID:$chain1:$chain2\t$id\t$alignmentLength\t$identities\t$interface\t";
               print FH "$asaChain1Single\t$asaChain1Complex\t$asaChain2Single\t$asaChain2Complex\n";
               close(FH);
            }
            else{
               open(FH,'>>',PISA_OUTPUT_HET) or die "Can not open/access '".PISA_OUTPUT_HET."'\n$!";
               select((select(FH), $|=1)[0]);
               print FH "$pdbID:$chain1:$chain2\t$id\t$alignmentLength\t$identities\t$interface\t";
               print FH "$asaChain1Single\t$asaChain1Complex\t$asaChain2Single\t$asaChain2Complex\n";
               close(FH);
            }
         }

#      $pm->finish;
      }
#   $pm->wait_all_children;
   }
}

sub pdb(){
   my ($HRsuperseeded, $cfg) = @_;
   my $HRpdbaa = $pdbaa->parseDunbrackFasta($cfg->findvalue('/cfg/path/pdbaa'));
   my $rootDir = $cfg->findvalue('/cfg/path/splitFiles');

   my $outHom = $cfg->findvalue('/cfg/path/pre_cull_hom');
   my $outHet = $cfg->findvalue('/cfg/path/pre_cull_het');

   open(FH,'>',$outHom) or die "Can not open/access '".$outHom."'\n$!";
      print FH "#PDB ID | \%seq id | alignment length | indentities  |  ASA interface  | ASA single chain 1  |  ASA complex chain 1 | ASA single chain 2  |  ASA complex chain 2  \n";
   close(FH);
   open(FH,'>',$outHet) or die "Can not open/access '".$outHet."'\n$!";
      print FH "#PDB ID | \%seq id | alignment length | indentities  |  ASA interface  | ASA single chain 1  |  ASA complex chain 1 | ASA single chain 2  |  ASA complex chain 2  \n";
   close(FH);
   my $sDirs = $pdbaa->parseDir($rootDir,'\w\w');
   for my $sdir (@{$sDirs}){
     my $subdir = $rootDir . $sdir;
     my $dirs = $pdbaa->parseDir($subdir,'\w\w\w\w');
     @{$dirs} = sort(@{$dirs});
     for my $dir (@{$dirs}){
       my $path = $rootDir . "$sdir/$dir/";
       my $ARfiles = $pdbaa->parseDir($path,'\w{4}_\w\w\.rsa');
       for my $rsaFile (@{$ARfiles}){
#      $pm->start and next;
#         next if(exists $HRoldRun->{$rsaFile});
         $rsaFile =~ /^(\w{4})_(\w)(\w)/;
         my $pdbID = $1;
         my $chain1 = $2;
         my $chain2 = $3;
         next if(exists $HRsuperseeded->{$pdbID});

         my $Lcomplex = $path . $rsaFile;
         my $Lsingle1 = $path . $pdbID . "_$chain1.rsa";
         my $Lsingle2 = $path . $pdbID . "_$chain2.rsa";

         if(not -e $Lcomplex){warn "No\t$Lcomplex";next}
         if(not -e $Lsingle1){warn "No\t$Lsingle1";next}
         if(not -e $Lsingle2){warn "No\t$Lsingle2";next}

         my $asaChain1Complex = &parseNaccessRsa($Lcomplex,$chain1);
         if(! $asaChain1Complex){print STDERR "NACCESS: No score single $pdbID:$chain1\t$Lcomplex\n"; next;}

         my $asaChain2Complex = &parseNaccessRsa($Lcomplex,$chain2);
         if(! $asaChain2Complex){print STDERR "NACCESS: No score single $pdbID:$chain2\t$Lcomplex\n"; next;}

         my $asaChain1Single = &parseNaccessRsa($Lsingle1,$chain1);
         if(! $asaChain1Single){print STDERR "NACCESS: No score single $pdbID:$chain1\t$Lsingle2\n"; next;}

         my $asaChain2Single = &parseNaccessRsa($Lsingle2,$chain2);
         if(! $asaChain2Single){print STDERR "NACCESS: No score single $pdbID:$chain2\t$Lsingle2\n"; next;}

         my $diff1 = $asaChain1Single - $asaChain1Complex;
         my $diff2 = $asaChain2Single - $asaChain2Complex;

         my $interface = ($diff1 + $diff2) / 2 ;
         $interface = sprintf("%d",$interface);
         if($interface > 0){
           my $pdbIDuc = uc($pdbID);
           my ($Fhom,$identities, $alignmentLength, $id) = &classifyInteraction($HRpdbaa,$pdbIDuc,$chain1,$chain2);
           if($Fhom == -1){next}


           if($Fhom == 1){
             open(FH,'>>',$outHom) or die "Can not open/access '".$outHom."'\n$!";
             select((select(FH), $|=1)[0]);
             print FH "$pdbID:$chain1:$chain2\t$id\t$alignmentLength\t$identities\t$interface\t";
             print FH "$asaChain1Single\t$asaChain1Complex\t$asaChain2Single\t$asaChain2Complex\n";
             close(FH);
           }
           else{
             open(FH,'>>',$outHet) or die "Can not open/access '".$outHet."'\n$!";
             select((select(FH), $|=1)[0]);
             print FH "$pdbID:$chain1:$chain2\t$id\t$alignmentLength\t$identities\t$interface\t";
             print FH "$asaChain1Single\t$asaChain1Complex\t$asaChain2Single\t$asaChain2Complex\n";
             close(FH);
           }
         }
       }
#      $pm->finish;
      }
#   $pm->wait_all_children;
   }
}

sub parsePisaOutput(){
   my %result;
   my ($Lfile1) = PISA_OUTPUT_HOM;
   my ($Lfile2) = PISA_OUTPUT_HET;

   open(FH,'<',PISA_OUTPUT_HOM) or die "Can not open/access '".PISA_OUTPUT_HOM."'\n$!";
      while(my $line = <FH>){

         if($line =~ /^(\w{4}):(\w):(\w)/){
         my @line = split("\t",$line);
         my $pdbID = $1;
         my $chain1 = $2;
         my $chain2 = $3;
         my $rsaFile =  $pdbID . "_" . $chain1.$chain2 .'.rsa';
         $result{$rsaFile}++;;
         }
      }
   close(FH);
   open(FH,'<',PISA_OUTPUT_HET) or die "Can not open/access '".PISA_OUTPUT_HET."'\n$!";
      while(my $line = <FH>){

         if($line =~ /^(\w{4}):(\w):(\w)/){
         my @line = split("\t",$line);
         my $pdbID = $1;
         my $chain1 = $2;
         my $chain2 = $3;
         my $rsaFile =  $pdbID . "_" . $chain1.$chain2 .'.rsa';
         $result{$rsaFile}++;;
         }
      }
   close(FH);
return(\%result);

}

=begin comment
#PDB ID | PDB reference chain | PISA original chain | PISA new chain
10gs  A  A  A
=cut
sub parsePisa2PdbMapping($){
   my %result;
   my ($Lfile) = @_;

   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
      while(my $line = <FH>){
         next if($line =~ /^#/);
         chomp($line);
         my @line = split("\t",$line);
         my $pdbID = $line[0];
         my $ref = $line[1];
         my $org = $line[2];
         my $new = $line[3];
         $result{$pdbID}{$new} = $ref;

      }
   close(FH);
return(\%result);

}

sub parseNaccessRsa($$){
   my ($Lfile, $chain) = @_;
   my $args;
   my $HRrsa = $pdbaa->parseNaccessRsa($Lfile);
   my $asa = $HRrsa->{'CHAIN'}->{$chain};
   if(! $asa){ return ()}
   if( ($asa =~ /[^0-9\.]/) || (! $asa) ){
      return();
   }
   return($asa);
}

sub classifyInteraction($\%){
   my ($HRpdbaa,$pdbID,$chain1,$chain2,$HRrefChains) = @_;
   my ($id, $identities,$alignmentLength,$longer);
   my $flagHom = 0;

   my $seq1 = $HRpdbaa->{'sequence'}->{$pdbID.$chain1};#print "\n".$pdbID.$chain1."\n";die;
   my $seq2 = $HRpdbaa->{'sequence'}->{$pdbID.$chain2};
   if(!$seq1 && $HRrefChains){
      my $pdbIDlc = lc($pdbID);
      my $refChain1 = $HRrefChains->{$pdbIDlc}->{$chain1};
#      print "1: $pdbID\t$refChain1\n"; die;
      $seq1 = $HRpdbaa->{'sequence'}->{$pdbID.$refChain1};
      if($refChain1){
         $seq1 = $HRpdbaa->{'sequence'}->{$pdbID.$refChain1};
      }
      else{
         if (!$seq1){warn "$pdbID"."$chain1 not found in pdbaa, check if superseeded"; return(-1)}
      }
   }
   if(!$seq2 && $HRrefChains){
      my $pdbIDlc = lc($pdbID);
      my $refChain2 = $HRrefChains->{$pdbIDlc}->{$chain2};
#      print "2: $pdbID\t$refChain2\t$chain2\n";
#     print Dumper($HRrefChains->{$pdbIDlc});  die;
      if($refChain2){
         $seq2 = $HRpdbaa->{'sequence'}->{$pdbID.$refChain2};
      }
      else{
         if (!$seq2){warn "$pdbID"."$chain2 not found in pdbaa, check if superseeded"; return(-1)}
      }
   }
   if (!$seq1){warn "$pdbID"."$chain1 not found in pdbaa, check if superseeded"; return(-1)}
   if (!$seq2){warn "$pdbID"."$chain2 not found in pdbaa, check if superseeded"; return(-1)}

   open(F1,">seq1") or die "Can not open/access 'seq1'\n$!\n" ;
   open(F2,">seq2") or die "Can not open/access 'seq2'\n$!\n" ;
      print F1 $seq1;
      print F2 $seq2;
   close(F1);
   close(F2);

   my $args = WATER . ' seq1 seq2 -gapopen 10.0 -gapextend 0.5 -sprotein1 -sprotein2 -outfile temp.out';
   system ($args) == 0 or die "Could not execute '$args' : $?";

   if (length($seq1) > length($seq2)){ $longer = length($seq1)}
   else{$longer = length($seq2)}
#more than 95% of the chain have to be mapped
   my $percent = ($longer / 100) * 95;


   open(FH,"<temp.out") or die "Can not open/access 'temp.out'\n$!\n";
LOOP:
   while(my $line = <FH>){
# Identity:    99/453 (21.9%)
      if($line =~ /^# Identity:/){
         $line =~ /\((.*)%\)$/;
         $id = sprintf("%d",$1);
         my @line = split(/\s+/,$line);
         if($line[2] !~ /^(\d+)\/(\d+)$/){print STDERR "Wrong format in local alignment output: '$line[2]'</water>\n";die;}
         $identities = $1;
         $alignmentLength = $2;
         last LOOP;
      }
   }
   close(FH);
# 95% or more identity and more than 95% of the length of the longer protein chain have to be covered
   if(!$alignmentLength){print "length: $pdbID\t$chain1\n$seq1\n$pdbID\t$chain2\n$seq2\n"; }
   if(!$percent){print "percent: $pdbID\t$chain1\n$seq1\n$pdbID\t$chain2\n$seq2\n"; }
   if(!$id){print "id: $pdbID\t$chain1\n$seq1\n$pdbID\t$chain2\n$seq2\n"; }
   if(  ($alignmentLength >= $percent) && ($id >= 95)  ){$flagHom = 1;}
   return($flagHom,$identities, $alignmentLength, $id);
}
