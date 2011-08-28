#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Pod::Usage;
use Parallel::ForkManager;
use File::Compare;
use Carp;

use Pdbaa;
use Constants;

pod2usage if ($ARGV[0] && ($ARGV[0] eq '-h' or $ARGV[0] eq '--help'));

=begin comment
use constant Constants->PDB_XRAY => '/home/tmp/BISC/update/pdb/pdbXray';
use constant Constants->PDBAA => '/home/tmp/BISC/update/pdbaa'; #23Aug2010
use constant Constants->PDBAA_OLD => '/home/tmp/BISC/update/pdbaa.26Apr2010'; 
use constant Constants->PDB_GET  => '/projects/compbio/bin/pdb-get -pdbdir  '; 
use constant Constants->PDB_STRUCTURES => '/home/tmp/BISC/update/pdb/structures/'; 
use constant Constants->PDB_STRUCTURES2 => '/home/tmp/BISC/update/pdb/structures2/'; 

use constant Constants->UPDATE_STATS => '/home/tmp/BISC/update/output/Bisc/updateStats'; 
use constant Constants->NEW_PDBSTRUCTURES => '/home/tmp/BISC/update/output/Bisc/newPDBstructures'; 
use constant Constants->OLD_PDBSTRUCTURES => '/home/tmp/BISC/update/output/Bisc/oldPDBstructures'; 
use constant Constants->OBSELETE_PDBSTRUCTURES => '/home/tmp/BISC/update/output/Bisc/obseletePDBstructures'; 
=cut
our $pdbaa = new Pdbaa;
select((select(STDOUT), $|=1)[0]); 

&main();

sub main(){
# Method returns hash that split information in FASTA comment line.
# This allows us to separate between X-ray and non-Xray methods
   my $HRpdbaa    = $pdbaa->parseDunbrackFasta(Constants->PDBAA);
   my $HRpdbaaOld = $pdbaa->parseDunbrackFasta(Constants->PDBAA_OLD);
   my $HRids      = &getXray($HRpdbaa);
   my $HRidsOld   = &getXray($HRpdbaaOld);
   my $HRnewSet   = &newSet($HRids,$HRidsOld);

   &printStats($HRids,$HRidsOld,$HRnewSet);
   &download($HRnewSet->{'new'});
   &download($HRnewSet->{'kept'});

}

sub printStats(\%\%\%){
   my ($HRids,$HRidsOld,$HRnewSet) = @_;

   open(FH,'>',Constants->UPDATE_STATS) or die "Can not open/access '".Constants->UPDATE_STATS."'\n$!";
      print FH 'Old pdbaa: '.Constants->PDBAA_OLD ."\n"; 
      print FH 'New pdbaa: '.Constants->PDBAA ."\n"; 

      print FH 'Old amount X-ray: '   . keys (%{$HRidsOld}) . "\n";
      print FH 'New amount X-ray: '   . keys (%{$HRids}) . "\n";
      print FH 'New complexes: '      . keys (%{$HRnewSet->{'new'}}) . "\n"; 
      print FH 'Complexes kept: '     . keys (%{$HRnewSet->{'kept'}}) . "\n"; 
      print FH 'Obselete complexes: ' . keys (%{$HRnewSet->{'obselete'}}) . "\n"; 
   close(FH);
   open(FH,'>',Constants->NEW_PDB_STRUCTURES) or die "Can not open/access '".Constants->NEW_PDB_STRUCTURES ."'\n$!";
      foreach my $pdbID (keys %{$HRnewSet->{'new'}}){
         print FH "$pdbID\n"; 
      }
   close(FH);  
   open(FH,'>',Constants->OLD_PDB_STRUCTURES) or die "Can not open/access '".Constants->OLD_PDB_STRUCTURES."'\n$!";
      foreach my $pdbID (keys %{$HRnewSet->{'kept'}}){
         print FH "$pdbID\n"; 
      }
   
   close(FH);  
   open(FH,'>',Constants->OBSELETE_PDB_STRUCTURES) or die "Can not open/access '".Constants->OBSELETE_PDB_STRUCTURES."'\n$!";
      foreach my $pdbID (keys %{$HRnewSet->{'obselete'}}){
         print FH "$pdbID\n"; 
      }
   
   close(FH);  

}
sub newSet(\%\%){
   my ($HRids,$HRidsOld) = @_;
   my $HRresult;
   foreach my $pdbID (keys %{$HRids} ){
      if(not exists $HRidsOld->{$pdbID} ){$HRresult->{'new'}->{$pdbID}++}
      else{$HRresult->{'kept'}->{$pdbID}++}
   }
   foreach my $pdbID (keys %{$HRidsOld}){
      if(not exists $HRids->{$pdbID}){$HRresult->{'obselete'}->{$pdbID}++}
   }
   return($HRresult);
}

sub getXray(\%){
   my ($HRpdbaa) = @_;
   my $HRout;

   foreach my $pdbID (sort keys %{$HRpdbaa->{'experiment'}}){
     next  if ($HRpdbaa->{'experiment'}->{$pdbID} ne 'XRAY');
     $pdbID = lc($pdbID);
     $HRout->{$pdbID}++;
   }
   return($HRout);
}
sub printXray(\%){
   my ($HRxray) = @_;
   open(FH,'>',Constants->PDB_XRAY) or die "Can not open/access '".Constants->PDB_XRAY."'\n$!";
      print FH "#X-ray structures: ".keys (%{$HRxray})."\n";
      foreach my $key (sort keys %{$HRxray}){
         print FH "$key\n";       
      } 
   close(FH);
}
sub download(){
   my ($HRids) = @_;
   foreach my $pdbID (keys %{$HRids} ){
      my $Lpdb = Constants->PDB_STRUCTURES .  "$pdbID.pdb";
      next if( -e $Lpdb);
      my $args = Constants->PDB_GET . Constants->PDB_STRUCTURES . " $pdbID";
      system($args) == 0 or die "$args";
      $args = 'gunzip ' . Constants->PDB_STRUCTURES . "$pdbID.pdb.gz";
      system($args) == 0 or die "$args";
      print "Finished: $pdbID\n"; 
   }
}

#      $pm->start and next; 
#      $pm->finish;
#   }
#   $pm->wait_all_children;
=head1 NAME

createPreBisc2010.pl


=head1 SYNOPSIS

createPreBisc2010.pl

=head1 AUTHOR

Thomas Juettemann <juettemann@gmail.com>


=head1 ARGUMENTS

=over 8

=item none

=back


=head1 OPTIONS

=over 8

=item none

=back

=head1 DESCRIPTION

B<createPreBisc2010.pl> compares 2 versions of Dunbracks pdbaa. All
PDB structures determined by X-ray are extracted from both files.
If the entry exists in the older file, but not in the newer one,
it will be marked obselete. If it exists in the newer one, but is 
missing in the older one it is listed as new. Entries that appear
in both are classified as kept.

=head1 BUGS

=over 4

=item None known
 
=back

=head1 TODO

=over 4

=item Nothing

=back

=head1 UPDATES

=cut


