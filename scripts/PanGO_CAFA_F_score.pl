#! /usr/local/bin/perl

#####
#####
##
##
###
#####
#####

# get command-line arguments
use Getopt::Std;
getopts('o:O:t:n:c:i:p:g:e:vVh') || &usage();
&usage() if ($opt_h);         # -h for help
$outFile = $opt_o if ($opt_o);    # -o for output files
$inFile = $opt_i if ($opt_i);     # -i for annotation with rank file
$existing = $opt_t if ($opt_t);   # -t for existing annotation file.
$gene_painted = $opt_g if ($opt_g);       # -g for painted gene file
$new = $opt_n if ($opt_n);        # -n for new annotations
$go_parent = $opt_p if ($opt_p);  # -p for go parent file
$evidence_count_out = $opt_O if ($opt_O); # -O for evidence count output
$errFile = $opt_e if ($opt_e);    # -e for (e)rror file (redirect STDERR)
$verbose = 1 if ($opt_v);         # -v for (v)erbose (debug info to STDERR)
$verbose = 2 if ($opt_V);         # -V for (V)ery verbose (debug info STDERR)


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


###################################
# Working on GO parents
###################################

my %go_parents;
my %go_child;
my %go;
#my %go_children;
open (GP, $go_parent);
while (my $line=<GP>){
    chomp $line;
   # next unless ($line=~/GO\:0003886/);
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
    
###########################################
#  Parse existing annotations
###########################################

my %existing;
open (EX, $existing);
while (my $line=<EX>){
    chomp $line;
    my ($id, $go, $aspect, $evidence, $type) = split(/\t/, $line);
    # next unless ($id eq 'Q16878');
    if (exists $uniprot_gene{$id}){
        my $longId = $uniprot_gene{$id};
        $existing{$longId}{$go}=$aspect;
        
        if (exists $go_parents{$go}){
            foreach my $parent (keys %{$go_parents{$go}}){
                next if ($parent=~/^GO:0005515|^GO:0005488|GO:0003674|GO:0008150|GO:0005575/);
                $existing{$longId}{$parent}=$aspect;
            }
        }
    }
}
close (EX);


#############################################################################
# Working on new annotations. In CAFA, this is considered as "ground truth"
#############################################################################

my %new;
open (NEW, $new);
while (my $line=<NEW>){
    chomp $line;
    my ($id, $go, $aspect, $evidence, $type) = split(/\t/, $line);

    if (exists $uniprot_gene{$id}){
        my $longId = $uniprot_gene{$id};
        next unless (exists $genes{$longId});    # only the ones in the curatable gene set
        next if (exists $existing{$longId}{$go});   # exclude pre-existing experimental annotations.
        $new{$aspect}{$longId}{$go}=1;
        $new{'all'}{$longId}{$go}=1;
    }
}
close (NEW);


################################################################################################
# Working on annotation prediction file
# These files already removed the existing annotations, so they are the predicted annotations.
################################################################################################

my %predicted;
my %id_in_both;
open (FH, $inFile);
while (my $line=<FH>){
    chomp $line;
    my ($id, $go, $name, $aspect, $score, $rank) = split(/\t/, $line);
    
    next unless (exists $new{$aspect}{$id});  # Only analyze genes that also have new annotations. This is how CAFA script does.
    $id_in_both{$aspect}{$id}=1;
    for ($i = 1; $i <=10; $i++){
        if ($i >= $rank){
            $predicted{$i}{$aspect}{$id}{$go}=1;
        }
        
    }
}
close (FH);

foreach my $aspect (keys %new){
    foreach my $id (keys %{$new{$aspect}}){
        next if (exists $id_in_both{$aspect}{$id});
        delete $new{$aspect}{$id};     # remove genes in new annotations that don't have predictions.
    }
}

######################
# Compare
######################

my @types = ("all", "mf", "bp", "cc");
foreach my $i (sort {$a<=>$b} keys %predicted){
    #next unless ($i eq 10);
    my %common;
    my %iba_more;
    foreach my $aspect (keys %{$predicted{$i}}){
        foreach my $id (keys %{$predicted{$i}{$aspect}}){
           # next unless ($id =~/Q9UHV7/);
            foreach my $go (keys %{$predicted{$i}{$aspect}{$id}}){
                if (exists $new{$aspect}{$id}{$go}){
                    $common{$id}{$aspect}{$go}=1;
                }else{
                    $iba_more{$id}{$aspect}{$go}=1;
                }
            }
        }
    }

    my %new_more;
    foreach my $aspect (keys %new){
        foreach my $id (keys %{$new{$aspect}}){
            #next unless ($id eq 'P09651');
            next unless (exists $predicted{$i}{$aspect}{$id}); #only IBA genes
            foreach my $go (keys %{$new{$aspect}{$id}}){
                next if (exists $common{$id}{$aspect}{$go});
                $new_more{$id}{$aspect}{$go}=1;
            }
        }
    }

    foreach my $id (keys %common){
        foreach my $aspect (keys %{$common{$id}}){
            my $go = join("\;", keys %{$common{$id}{$aspect}});
            #print "$id\t$aspect\t$go\n";
        }
    }
    my %gene_count;
    my %precision;
    my %recall;

    open (OUT, ">$output_file");
    foreach my $gene (keys %genes){
        foreach my $aspect (@types){
            
            my $n_common = &getCounts(\%common, $gene, $aspect);
            my $n_pred = &getCounts(\%iba_more, $gene, $aspect);
            my $n_exp = &getCounts(\%new_more, $gene, $aspect);
            
            my $common = join("\;", sort keys %{$common{$gene}{$aspect}});
            my $pred = join("\;", sort keys %{$iba_more{$gene}{$aspect}});
            my $exp = join("\;", sort keys %{$new_more{$gene}{$aspect}});
            
            my $total_pred = $n_common + $n_pred;
           
            my $total_exp = keys (%{$new{$aspect}{$gene}});  # all exp annotation of the aspect for the gene
            
            if ($total_pred >0){
                $gene_count{$aspect}{'pred'}{$gene}=1;
                my $precision = $n_common/$total_pred;
            
                $precision{$aspect}{$gene}{$precision}=1;
                
            }
            
            if ($total_exp >0){
                $gene_count{$aspect}{'exp'}{$gene}=1;
                my $recall = $n_common/$total_exp;
#
                $recall{$aspect}{$gene}{$recall}=1;
            }
        }
    }
    close (OUT);
    
    foreach my $type (@types){
        my $pred_count = keys (%{$gene_count{$type}{'pred'}});
        my $exp_count = keys (%{$gene_count{$type}{'exp'}});
        my $precision_total = &getTotal(\%precision, $type);
        my $recall_total = &getTotal(\%recall,$type);
        
        my $precision_ave = $precision_total/$pred_count;
        my $recall_ave = $recall_total/$exp_count;
        
        my $f_score;
        if (($precision_ave+$recall_ave) ==0){
            $f_score = 'na';
        }else{
            $f_score = (2*$precision_ave*$recall_ave)/($precision_ave+$recall_ave);
        }
        
        print "$i\t$type\t$precision_total\t$pred_count\t$precision_ave\t";
        print "$recall_total\t$exp_count\t$recall_ave\t";
        print "$f_score\n";
    }
}


##############################
# Subroutines
##############################

sub getCounts{
    my ($href, $gene, $aspect) = @_;
    
    my $n;
    if (exists $href->{$gene}->{$aspect}){
        $n = keys (%{$href->{$gene}->{$aspect}});
    }else{
        $n = 0;
    }
    return $n;
}

sub getTotal{
    my ($href, $t) = @_;
    
    my $n;
    foreach my $id (keys %{$href->{$t}}){
        foreach my $a (keys %{$href->{$t}->{$id}}){
            $n += $a;
            #print "$a\t$n\n";
        }
    }
    return $n;
}
