#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Pod::Usage;
use lib	'/home/tmp/BISC/libs';
use lib '/home/tmp/BISC/libs/lib/perl5/site_perl/5.8.8';
use Pdbaa;
use Parallel::ForkManager;

pod2usage if ($ARGV[0] && ($ARGV[0] eq '-h' or $ARGV[0] eq '--help'));

our $pdbaa = new Pdbaa;
our $pm = new Parallel::ForkManager(4);

use constant STARTDIR   => '/home/tmp/BISC/update/subcomplexes/pisa/';
#use constant STARTDIR  => '/home/tmp/BISC/update/subcomplexes/pdb/';
use constant NACCESS    => '/home/tmp/BISC/naccess2.1.1/naccess';
use constant LOGFILE    => '/home/tmp/BISC/scripts/2010/runNaccess.log';


&main();

sub main(){
   open(LOG,'>>',LOGFILE) or die "Can not open/access '".LOGFILE."'\n$!";
   select((select(LOG), $|=1)[0]);

   my $AR = $pdbaa->parseDir(STARTDIR,'\w\w');

   @{$AR} = sort(@{$AR});

   for my $sdir (@{$AR}){
      $pm->start and next;
      print LOG "STARTED: $sdir\n";
      my $path = STARTDIR . "$sdir/";
      my $ARfiles1 = $pdbaa->parseDir($path,'\w{4}_\w\.pdb');
      my $ARfiles2 = $pdbaa->parseDir($path,'\w{4}_\w\w\.pdb');

      for my $pdbFile (@{$ARfiles1}){
         my $Lpdb = $path . $pdbFile;
         my $LnaccessLog = $Lpdb;
         my $LnaccessRsa = $Lpdb;

         $LnaccessLog =~ s/pdb$/log/;
         $LnaccessRsa =~ s/pdb$/rsa/;

         if ( (not -T $LnaccessLog)  ||  (not -T $LnaccessRsa) ){
            chdir($path);
            &executeNaccess($Lpdb);
         }
      }
      for my $pdbFile (@{$ARfiles2}){
         my $Lpdb = $path . $pdbFile;
         my $LnaccessLog = $Lpdb;
         my $LnaccessRsa = $Lpdb;

         $LnaccessLog =~ s/pdb$/log/;
         $LnaccessRsa =~ s/pdb$/rsa/;

         if ( (not -T $LnaccessLog)  ||  (not -T $LnaccessRsa) ){
            chdir($path);
            &executeNaccess($Lpdb);
         }
      }
      print LOG "Finished $sdir\n";
      $pm->finish;
   }
   $pm->wait_all_children;
  close(LOG);
}
# execute NACCESS
sub executeNaccess($){
   my $Lfile = shift @_;
   my $args = NACCESS . " $Lfile";
   system($args) == 0 or die "Could not $args : $?";
   return()
}
=head1 NAME

runNAccess.pl


=head1 SYNOPSIS

runNAccess.pl

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

B<runNAccess.pl> executes NACCESS on all PISA and PDB subcomplexes. The directory of the subcomplexes is set as constant in the beginning of the file.

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


