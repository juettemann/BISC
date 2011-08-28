#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Pod::Usage;
use Pdbaa;
use Sql;

our $sql = new Sql;

pod2usage if ($ARGV[0] && ($ARGV[0] eq '-h' or $ARGV[0] eq '--help'));

use constant PDB_DIR => '/home/tmp/BISC/update/pdb/structures/';
use constant PISA_ORIGINAL => '/home/tmp/BISC/update/pisa/structures/original/';
use constant PISA_MODIFIED => '/home/tmp/BISC/update/pisa/structures/modified/';
use constant PDBAA => '/home/tmp/BISC/update/pdbaa.new';
use constant PISA_PDB_LOG => '/home/tmp/BISC/scripts/2010/test.pisa2pdb.log'; 
use constant PISA_PDB_ERR => '/home/tmp/BISC/scripts/2010/test.pisa2pdb.err'; 
use constant PISA_PDB_CHAINS => '/home/tmp/BISC/scripts/2010/test.pdb2pisaChains.log'; 
#use constant PISA_PDB_LOG => '/home/tmp/BISC/scripts/2010/pisa2pdb.log'; 
#use constant PISA_PDB_ERR => '/home/tmp/BISC/scripts/2010/pisa2pdb.err'; 
#use constant PISA_PDB_CHAINS => '/home/tmp/BISC/scripts/2010/pdb2pisaChains.log'; 


our $pdbaa = new Pdbaa;

&main();


#replaces chainId with new Chain
sub openLogs(){
   open(LOG,'>',PISA_PDB_LOG) or die "Can not open/access " . PISA_PDB_LOG  . "\n$!";
   open(LOG2,'>',PISA_PDB_CHAINS) or die "Can not open/access " . PISA_PDB_CHAINS  . "\n$!";
   open(ERRORLOG,'>',PISA_PDB_ERR) or die "Can not open/access " . PISA_PDB_ERR  . "\n$!";
   select((select(LOG), $|=1)[0]);
   select((select(LOG2), $|=1)[0]);
   select((select(ERRORLOG), $|=1)[0]);
   select((select(STDOUT), $|=1)[0]);
   print LOG "#PDB ID | PDB protein chains | PISA protein chains\n"; 
   print LOG2 "#PDB ID | PDB reference chain | PISA original chain | PISA new chain\n"; 
}
sub closeLogs(){

   close(ERRORLOG); 
   close(LOG); 
   close(LOG2); 
}
sub main(){
   &openLogs();
   my $dbh;
   $sql->connectBisc(\$dbh);
# Open logfiles and make them hot


   my $sthLog2 = $dbh->prepare(q{
         insert into pdb2pisaNew (pdbID, pdbRef, pisaOriginal, pisaNew )  values (?,?,?,?) ;
          }) or die $dbh->errstr;


#Files that contain the most protein chains in PISA
   my $ARpisaFiles = $pdbaa->parseDir(PISA_ORIGINAL,'pdb$');
   my $HRpdbaa = $pdbaa->parseDunbrackFasta(PDBAA);
#/home/tmp/BISC/update/pisa/structures/original/10gs.pisa.1.pdb
#   FILE:
   for my $pisaFile(@{$ARpisaFiles}){
#print STDOUT "$pisaFile\n"; 
      $pisaFile =~ /(\w{4})\.pisa\.(\d+)\.pdb$/;
#     next unless $pisaFile eq '/home/tmp/BISC/update/pisa/structures/original/1mam.pisa.1.pdb';
      my $pdbID = $1;
next unless($pdbID eq '2eh1');#
      my $serial = $2;
      my $pdbFile = PDB_DIR . "$pdbID.pdb";
      if(not -e $pdbFile){warn "$pdbFile does not exist";}
      my $Lout = PISA_MODIFIED . "$pdbID.pisa.modified.$serial.pdb";
      my $ucPdbID = uc($pdbID);
#print "Pisa file: $pisaFile\n"; 
#print "PDB  file: $pdbFile\n"; 
#print "Out: $Lout\n"; 

      my %result;
#get ATOM records and information about how they are mapped to the chains
      my ($HRpisa) = $pdbaa->parseATOMseq(PISA_ORIGINAL.$pisaFile);#print Dumper($HRpisa);die; 
      my ($HRpdb) = $pdbaa->parseATOMseq($pdbFile);
#determine amount of chains in PDB and PISA assemblies
      my $chainsPisa = keys(%{$HRpisa->{normal}});
      my $chainsPdb  = keys(%{$HRpdb->{normal}});
#some viruses have a huge amount of chains, exclude those
      if($chainsPisa > 62){
         print LOG "$pdbID\t$chainsPdb\t-1\n"; #
            next;
      }
#      print "$chainsPisa\t$chainsPdb\n"; 
#If we hit a chain that exists only in PISA, it needs a new chain ID that is not used for another chain in the PDB assembly
      my @range1 = ('A'..'Z');
      my @range2 = ('a'..'z');
      my @range3 = ('0'..'9');
      my @range = (@range1,@range2,@range3);

#      @range = @range[0 .. $chainsPisa];
      my %ids = map { $_ => 1 } @range;
#Remove the  IDs that exist in PDB to avoid identical chain identifiers
#      foreach my $chain (keys %{$HRpdb->{normal}}){
#         delete $ids{$chain} if(exists $ids{$chain});
#      } 
#sorting, 1st alphabetically, then numerical
      @range = &redoMerged(\%ids);
#print join(' ',@range),"\n"; 
#Identify all identical chains in PISA
      my $c = 0;
PISA1:
      foreach my $pisaXYZ (keys %{$HRpisa->{xyzToChains} }){
         print "$HRpisa->{xyzToChains}->{$pisaXYZ}\n"; 
         print "Her\n" if(exists $HRpdb->{xyzToChains}->{$pisaXYZ}); 
         next;
         my $pisaChain = $HRpisa->{xyzToChains}->{$pisaXYZ};
#print "PisaChain: $pisaChain\n"; 
#print "\nStart: $pisaChain\n"; 
         my $pisaSeq = $HRpisa->{normal}->{$pisaChain}->{sequence};
         my $refChain = $HRpdb->{referenceChain}->{$pisaSeq};
print "Ref Chain: $refChain\n"; 
         if(! $refChain){
            print ERRORLOG "No reference chain found for $pdbID\t$pisaChain\t$pisaFile\n"; 
            next;
         }
         if(exists $HRpdb->{xyzToChains}->{$pisaXYZ}){
            print $HRpdb->{xyzToChains}->{$pisaXYZ},$HRpisa->{xyzToChains}->{$pisaXYZ},"\n"; 
            print Dumper($HRpdb->{xyzToChains}->{$pisaXYZ}); 
            print Dumper($HRpisa->{xyzToChains}->{$pisaXYZ}); 
            if( ($HRpisa->{xyzToChains}->{$pisaXYZ} eq $HRpdb->{xyzToChains}->{$pisaXYZ}) ){
               $result{$pisaChain}{records} = $HRpisa->{normal}->{$pisaChain}->{records};
#print "Same Coordinates: $pisaChain\n"; 
               delete($ids{$pisaChain});
               @range = &redoMerged(\%ids);
               print LOG2 "$pdbID\t$refChain\t$pisaChain\t$pisaChain\n"; 
#               $sthLog2->execute($pdbID,$refChain,$pisaChain,$pisaChain) or die  $dbh->errstr;
               next PISA1;
            }
            else{
               my $newPisaChain = shift(@range);
               if(!$newPisaChain){
                  print LOG "$pdbID\t$chainsPdb\t-2\n"; #
                  next PISA1;
               }
               my $newRecord = &replaceChain($HRpisa->{normal}->{$pisaChain}->{records},$newPisaChain);
#print "Pisa Chain: $pisaChain\n"; 
#print "NewChain: $newPisaChain\n"; 
               $result{$newPisaChain}{records} = $newRecord;
               print LOG2 "$pdbID\t$refChain\t$pisaChain\t$newPisaChain\n"; 
               next PISA1;
            }
         }
         
      }
   
die;
=begin comment
print "PISA $chainsPisa\n"; 
foreach my $key (keys %{$HRpisa->{normal}}){
   print "$key "; 
} 
print "\n"; 
print "PDB $chainsPdb\n"; 
foreach my $key (keys %{$HRpdb->{normal}}){
   print "$key "; 
} 
print "\n"; 
=cut

#die;
#foreach my $key (sort keys %{$HRpisa->{normal}}){ print "$key\n"; } die;

#      print Dumper(\@range); die;
#foreach my $key (sort keys %ids){ print "$key "; } print "\n"; 

=begin comment
         Iterate through the PISA file
         Warn if sequence does not exist in PDB
         Assign PDB chain ID to PISA assembly
         Problem with homodimers. might not be able to get the same ID
=cut
   PISA:
      foreach my $pisaXYZ (keys %{$HRpisa->{xyzToChains} }){
#      print Dumper($HRpisa->{normal}->{$pisaChain}->{xyz}); 
         my $pisaChain = $HRpisa->{xyzToChains}->{$pisaXYZ};
print "PisaChain: $pisaChain\n"; 
#print "\nStart: $pisaChain\n"; 
         my $pisaSeq = $HRpisa->{normal}->{$pisaChain}->{sequence};
         my $refChain = $HRpdb->{referenceChain}->{$pisaSeq};
print "Ref Chain: $refChain\n"; 
         if(! $refChain){
            print ERRORLOG "No reference chain found for $pdbID\t$pisaChain\t$pisaFile\n"; 
            next;
         }
#print "RefChain: $refChain\n"; 
#print "PDBID: $ucPdbID\n"; 
#print "RefChain: $refChain\n"; 
#skip if Dunbrack does not recognise this as protein chain
         next PISA unless (exists $HRpdbaa->{'ids'}->{$ucPdbID}->{$refChain});
#print "Proper chain\n"; 
#Same coordinates present in PDB? 

# No similar coordinates found, try sequence comparison now
         my $pdbSeq = $HRpdb->{normal}->{$refChain}->{sequence};
         die if(!$pdbSeq);
# Sequence found. Now see if we have a matching chainID         
         if( exists $HRpdb->{seqToChains}->{$pisaSeq} ){
            my $size = scalar(@{$HRpdb->{seqToChains}->{$pisaSeq}}) ;
            for (my  $i = 0 ; $i < $size ; $i++){
               if(exists $HRpdb->{seqToChains}->{$pisaSeq}->[$i]){
                  if ($HRpdb->{seqToChains}->{$pisaSeq}->[$i] eq $pisaChain){
#print "Seq: $pisaChain\n"; 
                     $result{$pisaChain}{records} = $HRpisa->{normal}->{$pisaChain}->{records};
                     delete($HRpdb->{seqToChains}->{$pisaSeq}->[$i]);
                     @range = &redoMerged(\%ids);
                     print LOG2 "$pdbID\t$refChain\t$pisaChain\t$pisaChain\n"; 
#                  $sthLog2->execute($pdbID,$refChain,$pisaChain,$pisaChain) or die  $dbh->errstr;
                     next PISA;
                  }
               }
            }
            my $newPisaChain = shift(@range);
            if (!$newPisaChain ){
               
               my @newchains = sort keys %ids;
               $newPisaChain = $newchains[0];
               if (not defined $newPisaChain ){
                  print LOG "$pdbID\t$chainsPdb\t-1\n"; #
                     next PISA;
                  print "$pisaChain\n"; 
                  print "$pdbFile\n"; 
                  print Dumper($HRpdb); 
                  print Dumper(\@newchains); 
                  print Dumper(\%ids);
                  die; 
               }
               delete($ids{$newPisaChain});
               @range = &redoMerged(\%ids);
            }
            my $newRecord = &replaceChain($HRpisa->{normal}->{$pisaChain}->{records},$newPisaChain);
            $result{$newPisaChain}{records} = $newRecord;
#print "New Pisa: $newPisaChain\n"; 
            delete($ids{$newPisaChain});
            @range = &redoMerged(\%ids);
#PDB ID | PDB reference chain | PISA original chain | PISA new chain\n
            print LOG2 "$pdbID\t$refChain\t$pisaChain\t$newPisaChain\n"; 
#            $sthLog2->execute($pdbID,$refChain,$pisaChain,$newPisaChain) or die  $dbh->errstr;
            next PISA;

         }

# Same coordinates of first residue
#Check if other chain begins with those coordinates (ChainIDs just swapped in PISA
#no matching chain found, assign new one

      }
      

#PDB ID | PDB reference chain | PISA original chain | PISA new chain\n";
         
      my $proteinChainsPISA = keys %result;
      print LOG "$pdbID\t$chainsPdb\t$proteinChainsPISA\n"; 
#print "$Lout\n"; 
      open(FH,'>',$Lout) or die "Can not open/access '".$Lout."'\n$!";
      my $Frecords = 0;
         foreach my $chain (sort keys %result){
            if(defined $result{$chain}{records}){
               print FH "$result{$chain}{records}TER\n";
               $Frecords = 1;
            }
         }
         if($Frecords == 1){
            print FH "END\n"; 
         }
         else{
            print FH "Could not map chains in PISA assembly onto PDB\n"; 
         }
     close(FH);  
   }
#     print Dumper(\%result); 
   $sql->disconnectBisc(\$dbh); 
   &closeLogs;
}
#
sub redoMerged(\%){
   my ($HRids) = @_;
   my @range;
   my (@tempNumbers,@tempChar);

   foreach my $key (keys %{$HRids}){
      if($key =~ /\d/){
         push(@tempNumbers,$key);
      }
      elsif($key =~ /[a-zA-Z]/){
         push(@tempChar,$key);
      }
      else{die $key}
   }
   @range = (sort(@tempChar),sort{$a<=>$b}(@tempNumbers) );
  return(@range); 
}
sub replaceChain($$){
   if(scalar @_ != 2){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 2 needed\n";die;}
   my ($records,$newChain) = @_;
   die if(! $records);
   die if(not defined $newChain);
#print "'$newChain'\n"; 
   my $result;
   my @lines = split(/\n/,$records);
   for my $line(@lines){
      $line =~ /^(.{21})/;
      if(! $1){print Dumper(\@_); }
      $line =~ s/^(.{21})\w/$1$newChain/;
      my $lineStart = $1;
      if (! $lineStart) {die $line}
      $result .= "$line\n";
   }
   return($result);
}
#Check how many chains are assigned for each
#pisa and pdb have the same amount of chains
         
=begin commetn         
