package PdbMeta;

use 5.10.0;
use strict;
use warnings;

use Log::Log4perl qw(get_logger);
use Data::Dumper;
use XML::LibXML;
use Carp;

#use feature qw(switch say);

use constant UNIT  => {
  MONOMERIC      => 1,
  MONMERIC       => 1,
  DIMERIC        => 2,
  DIMER          => 2,
  TRIMERIC       => 3,
  TRIMER         => 3,
  TETRAMERIC     => 4,
  PENTAMERIC     => 5,
  HEXAMERIC      => 6,
  HEPTAMERIC     => 7,
  OCTAMERIC      => 8,
  NONAMERIC      => 9,
  DECAMERIC      => 10,
  UNDECAMERIC    => 11,
  DODECAMERIC    => 12,
  TRIDECAMERIC   => 13,
  TETRADECAMERIC => 14,
  PENTADECAMERIC => 15,
  HEXADECAMERIC  => 16,
  HEXADECAMRIC   => 16,
  HEPTADECAMERIC => 17,
  OCTADECAMERIC  => 18,
  NONADECAMERIC  => 19,
  EICOSAMERIC    => 20,
  '22MERIC'      => 22,
};

sub parseRemark350 {
  my ($self, $file) = @_;

  confess "File location is missing" if(!$file);
  confess "'$file' does not exists/not readable" if(! -T $file);

  my ($id, $reading) = (0,0);
  my $result;

  open(my $fh,'<', $file) or die "Can not open/access '$file'\n$!";
    while(my $line = <$fh>){
      next unless ($line =~ /^REMARK/);

      if (!$reading && $line =~ /^REMARK 350 BIOMOLECULE:\s+(\d+)/ ){
        $reading = 1;
        $id = $1;
      }
      elsif($reading &&
          $line =~ /AUTHOR DETERMINED BIOLOGICAL UNIT:\s+(\S*)/) {
        # 2yhn has nothing assigned
        $result->{$id}->{author} = 0 if(!$1);
        $result->{$id}->{author} = _getUnitValue($1);

      }
      elsif($reading &&
          $line =~ /SOFTWARE DETERMINED QUATERNARY STRUCTURE:\s+(\S*)/) {
        confess $line if(!$1);
        $result->{$id}->{software} = _getUnitValue($1);
      }
      elsif($reading &&
          $line =~ /SOFTWARE USED:\s+(\w+)/) {
        $result->{$id}->{method} = $1;
      }
      elsif($reading &&
          $line =~ /TOTAL BURIED SURFACE AREA:\s+(\d+)/) {
        $result->{$id}->{buried} = $1;
      }
      elsif($reading &&
          $line =~ /SURFACE AREA OF THE COMPLEX:\s+(\w+)/) {
        $result->{$id}->{surface} = $1;
      }
      elsif($reading &&
          $line =~ /CHANGE IN SOLVENT FREE ENERGY:\s+(\S+) /) {
        $result->{$id}->{energy} = $1;
      }
      elsif($reading &&
          $line =~ /^REMARK 350\s+$/){
        $reading = 0;
      }
      elsif($reading &&
          $line =~ /^REMARK 350 BIOMOLECULE:\s+(\d+)/){
        $id = $1;
      }
      else {
        $line =~ /REMARK\s*(\d{1,4})/;
        if(!$result && $reading){
          $result->{$id}->{author} = $id;
        }
        my $remark = $1;
        last if($remark > 350);
      }

    }
  close($fh);
  return($result);
}
###############################################################################
#
###############################################################################
sub _getUnitValue {
  my ($unit) = @_;
  confess "Unit missing" if(!$unit);
  $unit =~ s/\s//g;
  $unit =~ s/;//g;

  my $value;
  if($unit =~ /^(\d\d+)-?MERIC$/){
    $value = $1;
  }
  elsif($unit =~ /^(\w+)$/) {
    $value = UNIT->{$1};
    die "'$value' not defined in lookup" if(!$value);
  }
  else{die "Something odd with $unit"}
  return($value);
}
###############################################################################
#                           GetNumberOfBiomolecules
###############################################################################
# REMARK 300 BIOMOLECULE: 1, 2
sub getNumberOfBiomolecules {
  my ($self, $file) = @_;
  confess "File location is missing" if(!$file);
  confess "'$file' does not exists/not readable" if(! -T $file);

  my @molecules;
  my $found = 0;

  open(my $fh,'<', $file) or confess "Can not open/access '$file'\n$!";
  while(my $line = <$fh>){
    next unless ($line =~ /^REMARK 300 BIOMOLECULE:(.*)/);
    $found = 1;
    my $amount = $1;
    $amount =~ s/\s//g;
# make sure the numbers are consecutive
    if($amount =~ /,/){
      @molecules = split(/,/, $amount);
    }
    elsif(length($amount) == 1) {
      push(@molecules, $amount);
    }
    else{
      warn "Something wrong, line:\n$line\n";
    }
    last;
  }
  close($fh);

# No remark Biomolecule 300 exist, create it manuallt by parsing remark 350
  if(!$found){
  open(my $fh,'<', $file) or die "Can not open/access '$file'\n$!";
    while(my $line = <$fh>){
      next unless ($line =~ /^REMARK/);
      if ($line =~ /^REMARK 350 BIOMOLECULE:(.*)/){
        my $amount = $1;
        $amount =~ s/\s//g;
        die if($amount !~ /^[0-9]*$/);
        push(@molecules, $amount);
      }
      else{
        $line =~ /REMARK\s*(\d{1,4})/;
        my $remark = $1;
        last if($remark > 350);
      }
    }
  close($fh);

  }
  return(\@molecules);
}



###############################################################################
#                                       New                                   #
###############################################################################
# Constructor
###############################################################################
sub new {
  my ($class) = @_;

  my $self = {};

  bless($self, $class);
  return($self);
}
###############################################################################

1;
