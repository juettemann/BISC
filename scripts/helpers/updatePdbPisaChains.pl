#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Pod::Usage;
use Pdbaa;
use Sql;

pod2usage if ($ARGV[0] && ($ARGV[0] eq '-h' or $ARGV[0] eq '--help'));

use constant PDBAA => '/home/tmp/BISC/update/pdbaa.april26';

our $pdbaa = new Pdbaa;
our $sql = new Sql;


&main();
=begin comment
#PDB ID | PDB protein chains | PISA protein chains
10gs  2  2
10mh  1  1
117e  2  2
=cut
sub main(){
   my ($dbh);

   my $HRpdbaa = $pdbaa->parseDunbrackFasta(PDBAA);
   $sql->connectBisc(\$dbh);
   my $sth = $dbh->prepare(q{
         update pdbPisaChains set proteinPdb = ? where pdbID = ?; 
          }) or die $dbh->errstr;
   foreach my $pdbID (keys %{$HRpdbaa->{ids}}){
      my $chains = keys (%{$HRpdbaa->{ids}->{$pdbID}});
         $sth->execute($chains,$pdbID) or die  $dbh->errstr;
   }


   $sth->finish();

   $sql->disconnectBisc(\$dbh); 

}
=head1 NAME

script.pl


=head1 SYNOPSIS

script.pl

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

B<script.pl> 

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


