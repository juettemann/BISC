package Pdbaa;

use 5.10.0;

use strict;
use warnings;
use Data::Dumper;
use Carp;

use feature qw(say);


sub replaceChain {
  my ($self, $file) = @_;
  # flags
  my ($f_model, $f_chain) = (0,0);
  # counter for models in PDB file
  my ($model, $endMdl)    = (0,0);

  # Read in file
  open(my $fh,'<',) or confess "Can not open/access '$file'\n$!";
  while(my $line = <$fh>){
    my @file = <$fh>;
  }
  close($fh);

  foreach my $line (@file){
    chomp($line);

    if ($line =~ /^MODEL/){
      $model++;
    }
    elsif ($line =~ /^ENDMDL/){
      $endMdl++;
    }
    elsif ($line =~ /^TER/){
    }
    elsif ( (($line =~ /^ATOM/) || ($line =~ /^HETATM/)) &&
        ){

      my $serial      = substr($line,7,4);  $serial      =~ s/^\s+|\s+$//g; #0
      my $atomName    = substr($line,13,3); $atomName    =~ s/^\s+|\s+$//g; #1
      my $altLocation = substr($line,16,1); $altLocation =~ s/^\s+|\s+$//g; #2
      my $residueName = substr($line,17,3); $residueName =~ s/^\s+|\s+$//g; #3
      my $chain       = substr($line,21,1); $chain       =~ s/^\s+|\s+$//g; #4
      my $resNo       = substr($line,23,3); $resNo       =~ s/^\s+|\s+$//g; #5
      my $x           = substr($line,31,7); $x           =~ s/^\s+|\s+$//g; #6
      my $y           = substr($line,39,7); $y           =~ s/^\s+|\s+$//g; #7
      my $z           = substr($line,47,7); $z           =~ s/^\s+|\s+$//g; #8
      my $occupancy   = substr($line,55,5); $occupancy   =~ s/^\s+|\s+$//g; #9
      my $tempFactor  = substr($line,61,5); $tempFactor  =~ s/^\s+|\s+$//g; #10
      my $element     = substr($line,77,1); $element     =~ s/^\s+|\s+$//g; #11
      my $charge      = substr($line,79,1); $charge      =~ s/^\s+|\s+$//g; #12
      $result{$endMdl}{$chain}  .= "$line\n";
    }
  }

}

sub new {
  my $class = shift;
  my $self = {};
  my %param = @_;
  bless($self, $class);
  return $self;
}
1;
