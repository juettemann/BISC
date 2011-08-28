#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Pod::Usage;
use Pdbaa;
use Sql;


pod2usage if ($ARGV[0] && ($ARGV[0] eq '-h' or $ARGV[0] eq '--help'));

use constant PDB2PISA => '/home/tmp/BISC/scripts/2010/pdb2pisaChains.log'; 


our $pdbaa = new Pdbaa;
our $sql = new Sql;

&main();
=begin comment
==> test.pdb2pisaChains.log <==
#PDB ID | PDB reference chain | PISA original chain | PISA new chain
2eh1  B  F  B
2eh1  A  A  A
=cut
=begin comment
mysql> desc pdb2pisaNew;
+--------------+---------+------+-----+---------+----------------+
| Field        | Type    | Null | Key | Default | Extra          |
+--------------+---------+------+-----+---------+----------------+
| id           | int(11) | NO   | PRI | NULL    | auto_increment | 
| pdbID        | char(4) | NO   |     |         |                | 
| pdbRef       | char(1) | NO   |     |         |                | 
| pisaOriginal | char(1) | NO   |     |         |                | 
| pisaNew      | char(1) | NO   |     |         |                | 
+--------------+---------+------+-----+---------+----------------+
=cut

sub main(){
   my $dbh;
   $sql->connectBiscTest(\$dbh);
   $dbh->do('truncate pdb2pisa');
   my $sth = $dbh->prepare(q{
         insert into pdb2pisa ( pdbID, pdbRef, pisaOriginal, pisaNew)  values (?,?,?,?) ;
          }) or die $dbh->errstr;
   my $HR;


   open(FH,'<',PDB2PISA) or die "Can not open/access '".PDB2PISA."'\n$!";
      while(my $line = <FH>){
         next if($line =~ /^#/);
         chomp($line);
         $line =~ /^(\w{4})\t(\w)\t(\S*)\t(\w)/;
         $sth->execute($1,$2,$3,$4) or die  $dbh->errstr;
      }
      close(FH);

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



