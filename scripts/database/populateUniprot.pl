#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Pod::Usage;
use DBI;
use Pdbaa;
use Sql;
#DBI->trace( 2 );

pod2usage if ($ARGV[0] && ($ARGV[0] eq '-h' or $ARGV[0] eq '--help'));

our $pdbaa = new Pdbaa;
our $sql   = new Sql;

use constant UNIPROT_ACCS => '/home/tmp/BISC/uniprot/modeltargets.accs';
use constant UNIPROT_FASTA => '/home/tmp/BISC/uniprot/modeltargets.fasta';
use constant UNIPROT_TXT => '/home/tmp/BISC/uniprot/modeltargets.txt';

&main();
=begin comment
In order to get smaller tables and quicker queries, only sequences in modeltargets are stored in this table
Read UniProt text file
Read UniProt fasta file
Get all ids from BISC-MI
update uniprot
=cut

sub main(){
   my $dbh;
   $sql->connectBisc(\$dbh);
   my $HRaccs = getAccs($dbh);
#   die;
#   print Dumper($HRaccs); 
   my $HRupTxt   = $pdbaa->parseUniprotText(UNIPROT_TXT);
   print "Got txt\n"; 
   my $HRupFasta = $pdbaa->parseUniprotFasta(UNIPROT_FASTA);

   print "Got Fasta\n"; 
   &populateUniprot($dbh,$HRaccs,$HRupTxt,$HRupFasta);
   $sql->disconnectBisc(\$dbh); 
   
}
sub populateUniprot(){
   if(scalar @_ != 4){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 4 needed\n";die;}
   my ($dbh,$HRaccs,$HRupTxt,$HRupFasta) =  @_;
#  print Dumper($HRentries);
   my ($lenAcc,$lenCline,$lenSeq,$lenTaxid) = (0,0,0,0);

   my $sth = $dbh->prepare(q{
      insert into uniprot (accession,taxid,sequence,fasta,description) values (?,?,?,?,?);
         }) or die $dbh->errstr;

   my $rows_affected = $dbh->do("truncate uniprot");
   foreach my $acc (sort keys %{$HRaccs}){
      my $taxid = $HRupTxt->{taxid}->{$acc};
      my $seq = $HRupFasta->{sequence}->{$acc};
      my $fasta = $HRupFasta->{fasta}->{$acc};
      my $description = $HRupTxt->{description}->{$acc};
#      if($acc eq '')
      if(!$taxid){ 
         warn "$acc\n"; 
         next;
#         die;
         #print Dumper($HRupTxt->{taxid});
      }
die $acc if(!$seq);
      $sth->execute($acc,$taxid,$seq,$fasta,$description) or die $dbh->errstr ;
      if(length($description) > $lenCline){$lenCline = length($description)}
      if(length($seq) > $lenSeq){$lenSeq = length($seq)}
      if(length($acc) > $lenAcc){$lenAcc = length($acc)}
      if(length($taxid) > $lenTaxid){$lenTaxid = length($taxid)}

   }
   print "Out\n"; 
   $sth->finish();

   print "Length\nAcc: $lenAcc\nCline:$lenCline\nSeq:$lenSeq\nTaxid:$lenTaxid\n\n"; 
   return();
}



sub getAccs($\%){
   if(scalar @_ != 1){print STDERR "Argument list length: " . scalar(@_) . " argument(s) passed, 1 needed\n";die;}
   my ($dbh) = @_;
   my $HRaccs;

   my $sthHom = $dbh->prepare(q{
         select accession1,accession2  from biscMiHom;
         }) or die $dbh->errstr;
   my $sthHet = $dbh->prepare(q{
         select accession1,accession2  from biscMiHet;
         }) or die $dbh->errstr;
   $sthHom->execute() or die  $dbh->errstr;
   
   while(my ($acc1,$acc2)  = $sthHom->fetchrow_array()){
      $HRaccs->{$acc1}++;
      $HRaccs->{$acc2}++;
   }
   $sthHom->finish();

   $sthHet->execute() or die  $dbh->errstr;
   while(my ($acc1,$acc2)  = $sthHet->fetchrow_array()){
      $HRaccs->{$acc1}++;
      $HRaccs->{$acc2}++;
   }
   $sthHet->finish();
   open(FH,'>',UNIPROT_ACCS) or die "Can not open/access '".UNIPROT_ACCS."'\n$!";
      foreach my $acc (keys %{$HRaccs}){
         print FH "$acc\n"; 
      }
   close(FH);  


   return($HRaccs);
}





=head1 NAME

populateUniprot.pl


=head1 SYNOPSIS

populateUniprot.pl

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

B<populateUniprot.pl> populates uniprot table in BISC MySQL database.


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


