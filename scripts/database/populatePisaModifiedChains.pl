#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Pdbaa;
use Sql;
use Parallel::ForkManager;


use constant ASSEMBLIES => '/home/tmp/BISC/update/pisa/lists/assemblies/';
use constant PDB_STRUCTURES => '/home/tmp/BISC/update/pdb/structures/';
use constant PISA_MODIFIED => '/home/tmp/BISC/update/pisa/structures/modified/';
use constant PDBAA  => '/home/tmp/BISC/update/pdbaa';
use constant PISA2PDB => '/home/tmp/BISC/scripts/2010/test.pisa2pdb.log'; 
use constant PISA2PDB_OUT => '/home/tmp/BISC/scripts/2010/pisa2pdb_DB'; 
use constant PISACFG1 => '/home/tmp/BISC/pisa/mmdb/pisa.bisc1.cfg'; 
use constant PISACFG2 => '/home/tmp/BISC/pisa/mmdb/pisa.bisc2.cfg';

our $sql = new Sql;
our $pdbaa = new Pdbaa;
our $pm = new Parallel::ForkManager(4);

&main();

sub main(){
   
   my $ARpisa = $pdbaa->parseDir(PISA_MODIFIED,'\w{4}\.pisa.modified\.\d*\.pdb');
   print Dumper($ARpisa); 
}
