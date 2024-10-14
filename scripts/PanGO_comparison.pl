#! /usr/local/bin/perl

#####
#####
##
##   -i preidcction_PanGO file
##   -a prediction_[method]
##   -g PanGO_curated_genes
##   -p go_all_parents file
###
#####
#####

# get command-line arguments
use Getopt::Std;
getopts('o:O:i:n:a:d:p:g:e:vVh') || &usage();
&usage() if ($opt_h);         # -h for help
$outFile = $opt_o if ($opt_o);    # -o for
$inFile = $opt_i if ($opt_i);     # -i for IBA file
$inDir = $opt_d if ($opt_d);      # -d for input Dir for GO annotation files
$annotation = $opt_a if ($opt_a); # -a for the specific annotation file for each method.
$gene_painted = $opt_g if ($opt_g);       # -g for gene.dat file
$go_parent = $opt_p if ($opt_p);   # -p for go parent file.
$errFile = $opt_e if ($opt_e);    # -e for (e)rror file (redirect STDERR)
$verbose = 1 if ($opt_v);         # -v for (v)erbose (debug info to STDERR)
$verbose = 2 if ($opt_V);         # -V for (V)ery verbose (debug info STDERR)

$| = 1;

############################################
# Parse the painted gene file
############################################
my %genes;
my %uniprot_gene;  # uniprotID to longID mapping
open (PG, $gene_painted);
while (my $line=<PG>){
    chomp $line;
    $genes{$line}=1;
    my $uniprot = (split(/\=/, $line))[2];
    $uniprot_gene{$uniprot}=$line;
}
close (PG);


####################################
# Working on all GO parent file.
####################################

my %go_parents;
my %go_child;
my %go;

open (GP, $go_parent);
while (my $line=<GP>){
    chomp $line;

    my ($child, $parent, $relation, $foo)=split(/\t/, $line);
    
    my ($child_name, $child_id);
    
    if ($child=~/^(.+)\((GO\:\d+)\)$/){
        $child_name=$1;
        $child_id = $2;
    }
    
    my ($parent_name, $parent_id);
    if ($parent=~/^(.+)\((GO\:\d+)\)$/){
        $parent_name = $1;
        $parent_id = $2;
    }
##
    my $aspect;
    if ($foo=~/^m/){
        $aspect = 'mf';
    }elsif ($foo=~/^b/){
        $aspect = 'bp';
    }elsif ($foo=~/^c/){
        $aspect = 'cc';
    }
    $go_parents{$child_id}{$parent_id}=1;
    $go{$child_id}="$aspect\t$child_name";
    $go{$parent_id}="$aspect\t$parent_name";
    $go_child{$parent_id}{$child_id}=1;
    # print "$line\n";
}
close (GP);


#########################################################################
# Working on prediction_PanGO file (existing annotations were removed already
#########################################################################

my %IBA;
open (IBA, $inFile);
while (my $line=<IBA>){
    chomp $line;
    
    my ($id, $go, $name, $aspect, $score, $rank, $type) = split(/\t/, $line);
    next unless ($type eq 's');  # This is the direct annotation.
    $IBA{$aspect}{$id}{$go}=1;
}
close (IBA);


################################################################################################
# Working annotation prediction file
# These files already removed the existing annotations, so they are the predicted annotations.
################################################################################################

my %predicted;
my %id_in_both;
open (FH, $annotation);
while (my $line=<FH>){
    chomp $line;
    my ($id, $go, $name, $aspect, $score, $rank) = split(/\t/, $line);
    next unless (exists $genes{$id});
    
    next unless (exists $IBA{$aspect}{$id});  # include genes in IBA only.
    $id_in_both{$aspect}{$id}=1;
    for ($i = 1; $i <=10; $i++){
        if ($i >= $rank){
            $predicted{$i}{$aspect}{$id}{$go}=1;
        }
    }
}
close (FH);

foreach my $aspect (keys %IBA){
    foreach my $id (keys %{$IBA{$aspect}}){
        next if (exists $id_in_both{$aspect}{$id});
        delete $IBA{$aspect}{$id};  # remove genes that are not in the prediction.
    }
}


my @aspects = ('mf', 'bp', 'cc');
print "Threshold pools\tOntology Aspect\tSame term\tMore specific term in IBA\tMore specific term in method\tMore general term in IBA\tMore general term in method\tOnly in IBA\tOnly in method\n";
foreach my $i (sort {$a<=>$b} keys %predicted){
    print STDERR "Working on pool $i\n";
    foreach my $aspect (@aspects){
        my %common;
        my %iba_parent;
        my %iba_child;
        my %iba_only;
        my %method_child;
        my %method_parent;
        my %method_only;
        foreach my $id (keys %{$predicted{$i}{$aspect}}){
            my %non_redundant_go;
            foreach my $go (keys %{$predicted{$i}{$aspect}{$id}}){
                my $check;
                foreach my $go1 (keys %{$predicted{$i}{$aspect}{$id}}){
                    if (exists $go_parents{$go1}{$go}){
                        $check++;
                        last;
                    }
                }
                next if ($check);
                $non_redundant_go{$go}=1;
            }
         
            if (exists $IBA{$aspect}{$id}){
                foreach my $go (keys %{$IBA{$aspect}{$id}}){
                    if (exists $non_redundant_go{$go}){
                        $common{"$id\t$go"}=1;
                    }else{
                        my $check1;
                        foreach my $go1 (keys %non_redundant_go){
                            if (exists $go_parents{$go}{$go1}){
                                $iba_child{"$id\t$go"}=1;
                                $method_parent{"$id\t$go1"}=1;
                                $check1++;
                            }elsif (exists $go_child{$go}{$go1}){
                                $iba_parent{"$id\t$go"}=1;
                                $method_child{"$id\t$go1"}=1;
                                $check1++
                            }
                        }
                        next if ($check1);
                        $iba_only{"$id\t$go"}=1;
                        
                    }
                }
            }
            
            foreach my $go (keys %non_redundant_go){
                next if (exists $common{"$id\t$go"});
                next if (exists $method_parent{"$id\t$go"});
                next if (exists $method_child{"$id\t$go"});
                $method_only{"$id\t$go"}=1;
            }
        }
        
        foreach my $id (keys %{$IBA{$aspect}}){
            foreach my $go (keys %{$IBA{$aspect}{$id}}){
                my $foo = "$id\t$go";
                next if (exists $common{$foo});
                next if (exists $iba_child{$foo});
                next if (exists $iba_parent{$foo});
                $iba_only{$foo}=1;
            }
        }
        
        my $total_common = keys (%common);
        my $total_iba_child = keys (%iba_child);
        my $total_method_child = keys (%method_child);
        my $total_iba_parent = keys (%iba_parent);
        my $total_method_parent = keys (%method_parent);
        my $total_iba_only  = keys (%iba_only);
        my $total_method_only = keys (%method_only);
        
       
        print STDERR "Printing $aspect of pool $i\n";
        print "$i\t$aspect\t$total_common\t$total_iba_child\t$total_method_child\t$total_iba_parent\t$total_method_parent\t$total_iba_only\t$total_method_only\n";
    }
}

