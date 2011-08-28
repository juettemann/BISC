#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Pod::Usage;
use Pdbaa;
use Sql;

pod2usage if ($ARGV[0] && ($ARGV[0] eq '-h' or $ARGV[0] eq '--help'));

our $sql = new Sql;
our $pdbaa = new Pdbaa;




&main();

sub main(){
   $sql->connectBisc(\$dbh);
   my $sth = $dbh->prepare(q{
         update pdbPisaChains set proteinPisa = ? where pdbID = ?; 
          }) or die $dbh->errstr;

   open(FH,'<',$Lfile) or die "Can not open/access '".$Lfile."'\n$!";
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


