#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Pdbaa;

use constant PDBAA  => '/home/tmp/BISC/update/pdbaa.new';
use constant LATEST => '/home/tmp/BISC/update/s2c/Latest/';
#use constant LATEST => '/home/tmp/Latest/';
use constant OUTPUT => '/home/tmp/BISC/update/pdbaa_atom.test';
#use constant OUTPUT => '/home/s0571283/3eb7.atom.fasta';

our $pdbaa = new Pdbaa;
our $HRpdbaa = $pdbaa->parseDunbrackFasta(PDBAA);
#print Dumper($HRpdbaa->{'commentLine'});die; 

my $HRatomsequence = &createAtomSequence();
&printOutput($HRatomsequence);

sub printOutput(){
   my $HRatomsequence = shift @_;
   open(FH,'>',OUTPUT) or die "Can not open/access '".OUTPUT."'\n$!"; 
      foreach my $pdbAndChainID (sort keys %{$HRatomsequence}){
         my $header;
         if(not exists $HRpdbaa->{'commentLine'}->{$pdbAndChainID}){
            print STDERR "No header for $pdbAndChainID, either DNA, obsolete, to short or not using most recent pdbaa etc.\n";
            $header = '>' . $pdbAndChainID;
         }
         else{ $header = $HRpdbaa->{'commentLine'}->{$pdbAndChainID};}
         my $sequence = $HRatomsequence->{$pdbAndChainID};
            $sequence =~ s/(.{80}|.{1,79}$)/$1\n/g;
         print FH "$header\n";
         print FH "$sequence\n";
      }
   close(FH);
#%print Dumper($HRatomsequence); die;
}

sub createAtomSequence($$){
   opendir(DIR,LATEST);
      my @files = grep(/\.sc$/,readdir(DIR));
   closedir(DIR);
#print Dumper(\@files);die;
my %stdRes = ( 'ALA' , '1' ,'ARG' , '1' ,'ASN' , '1' ,'ASP' , '1' ,'CYS' , '1' ,'GLU' , '1' ,'GLN' , '1' ,'GLY' , '1' ,'HIS' , '1' ,'ILE' , '1' ,'LEU' , '1' ,'LYS' , '1' ,'MET' , '1' ,'PHE' , '1' ,'PRO' , '1' ,'SER' , '1' ,'THR' , '1' ,'TRP' , '1' ,'TYR' , '1' ,'VAL' , '1','MSE','1' );


   my (%atomSequence);
   foreach my $file (@files){
      my $pdbID = substr($file,0,4);
      my $filename = LATEST . $file;
      open(SC,"<$filename") or die "Cannot open/access '".$filename."'\n$!";
         while(my $line = <SC>){
            if($line =~ /^SEQCRD/){
               my $chain   = substr($line,7,1);
               my $residue = substr($line,9,1);
               my $name    = substr($line,15,3);
               if(not exists $stdRes{$name} ){$residue = 'X';}
               if( $name   =~ /---/){$residue = 'X';}
               $pdbID =~ tr/[a-z]/[A-Z]/;
               $atomSequence{$pdbID.$chain} .= $residue;
            }
         }   
      close(SC);
   }
   return(\%atomSequence);
}
