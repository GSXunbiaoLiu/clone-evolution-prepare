#!/usr/bin/perl -w
use strict;
use lib '/GPFS01/home/zhouyh/perl5/lib/perl5/';
use FindBin '$Bin';
use Getopt::Long;
use DateTime;
use File::Basename;
use Cwd;
sub usage
{
        print STDERR <<USAGE;
==============================================================================
Description    
Options
        -outdir <s>: pathway of outdir
        -info <s>: info file
        -bamdir <s>: pathway of bam files
        -segdir <s>: pathway of seg files
        -fa : pathway of reference file
	-mafdir : pathway of maf files
	-ex : list of infos to exclude
        -h|?|help : Show this help
==============================================================================
USAGE
}

my ($help,$cpu,$nodes,$info,$mafdir,$bamdir,$segdir,$outdir,$fa,$ex);
GetOptions(
        "h|?|help"=>\$help,
        "info=s"=>\$info,
        "mafdir=s"=>\$mafdir,
        "fa=s"=>\$fa,
        "segdir=s"=>\$segdir,
        "outdir=s"=>\$outdir,
	"bamdir=s"=>\$bamdir,
        "cpu=s"=>\$cpu,
        "nodes=s"=>\$nodes,
	"ex=s"=>\$ex
);
if( defined($help) || !defined ($info)|| !defined ($mafdir)|| !defined ($segdir)){
        &usage;
        exit 0;
}

open (INF,"$info");
$fa ||= "/GPFS01/home/njsh/databases/GSCAP/hg19.fa";
my %info;
$outdir ||= getcwd();
$cpu = 5;
$nodes = "all";
`mkdir $outdir/absolute_result` unless (-e "$outdir/absolute_result");
`mkdir $outdir/absolute_result/purity` unless (-e "$outdir/absolute_result/purity");
open (ABS,">$outdir/absolute_result/absolute.sh");
while(<INF>){
	chomp;
	my ($patient,$name) = split /\,/,$_;
	push @{$info{$patient}},$name;
	my $absolute_sh = absolute ($name,$segdir,$mafdir,"$outdir/absolute_result","$outdir/absolute_result/purity");
	print ABS "$absolute_sh\n";
}
close INF;
close ABS;

chdir ("$outdir/absolute_result");
`bsubjobs absolute.sh`;

chdir ("$outdir");
`mkdir $outdir/combine_maf` unless (-e "$outdir/combine_maf");
`mkdir $outdir/GetBasecounts` unless (-e "$outdir/GetBasecounts");
`mkdir $outdir/PyClone_input` unless (-e "$outdir/PyClone_input");
foreach my $patient (keys %info){
	open (PATIENT,">$outdir/combine_maf/$patient.somatic.filter.maf.txt");
	my $combined_maf = cat_maf ($info{$patient},$mafdir);
	print PATIENT "$combined_maf";
	close PATIENT;
        open (ADD,">$outdir/combine_maf/$patient.somatic.filter.maf_add_chr.txt");
	my $add_chr_maf = add_chr ("$outdir/combine_maf/$patient.somatic.filter.maf.txt");
	print ADD "$add_chr_maf";
	close ADD;
	open (GBC,">$outdir/GetBasecounts/$patient.getbasecount.sh");
	my $GetBasecounts = GetBasecounts ($info{$patient},$bamdir,$fa,$patient,"$outdir/combine_maf/$patient.somatic.filter.maf_add_chr.txt","$outdir/GetBasecounts");
	print GBC "$GetBasecounts\n";
	close GBC;
	chdir ("$outdir/GetBasecounts");
	&pid_stat ("$patient.getbasecount.sh");
	open (ADD_HGVS,">$patient.GetBaseCountsMultiSample.txt");
	my $add_HGVS = add_HGVSp_Short("$outdir/combine_maf/$patient.somatic.filter.maf_add_chr.txt","$outdir/GetBasecounts/$patient.GetBaseCountsMultiSample.omaf.txt");
	print ADD_HGVS "$add_HGVS\n";
	close ADD_HGVS;
	chdir ("$outdir");
	Pyclone_format ("$outdir/GetBasecounts/$patient.GetBaseCountsMultiSample.txt",$segdir,"$outdir/PyClone_input");
}
if ($ex){
	`mkdir $outdir/PyClone_filter_input`;
	my @infile = `ls $outdir/PyClone_input/*`;
	foreach my $infile (@infile){
		chomp $infile;
		my $name = basename ($infile);
	        my @name = split /\./,$name;
	        my $file_name = $name[0];
#		print "$file_name\n";
		open (OUT,">$outdir/PyClone_filter_input/$file_name.PyClone.filter.input");
		my $filter_out = exclude_list ($infile,$ex);
		print OUT "$filter_out";
		close OUT;
	}
}

sub exclude_list {
	my ($file,$exclude) = @_;
	my $out_exclude;
	my %exclude;
	open (IN,"$exclude");
	while (<IN>){
		chomp;
		my @all = split /\t+/,$_;
		my $vinf;
		if ($all[-2] ne "."){
			$vinf = "$all[0]--$all[-2]";
		} else {
			$vinf = "$all[0]--$all[-3]";
		}
		push @{$exclude{$all[-1]}},$vinf;
	}
	close IN;
	my $name = basename ($file);
	my @name = split /\./,$name;
	my $file_name = $name[0];
	if (exists $exclude{$file_name}){
		my @values = @{$exclude{$file_name}};
		open (IN,"$file");
		my $head = <IN>;
		$out_exclude .= $head;
		while (<IN>){
			my $i = 0;
			foreach my $value (@values){
				if ($_ =~ /$value/){
					$i = 1;
					next;
				}
			}
			if($i == 0){
				$out_exclude .= $_;
			}
		}
	 }else{
		open (IN,"$file");
		while (<IN>){
			$out_exclude .= $_;
		}
	}
	close IN;
	return $out_exclude;
}

sub add_HGVSp_Short {
	my ($combine_maf,$omaf) = @_;
	my %mut;
	open (IN,"$combine_maf");
	my %title;
        my $head = <IN>;
        if($head =~ /^#/){
                $head = <IN>;
        }
#        print "$head";
        die "ERROR:combine_maf file must have title!\n" if ($head !~ /Hugo_Symbol/);
        chomp $head;
        my @title = split /\t/,$head;
        for my $i (0..$#title){
                $title{$title[$i]} = $i;
        }
	while (<IN>){
		my @all = split /\t/,$_;
			my $k = "$all[$title{Chromosome}]\t$all[$title{Start_Position}]\t$all[$title{Variant_Classification}]";
#			print "$k\n";
			my $v = "$all[$title{HGVSc}]\t$all[$title{HGVSp_Short}]\t$all[$title{HGVSc}]";
			$mut{$k} = $v;
	}
	close IN;
	open (IN,"$omaf");
	my $add_HGVSp_Short_out;
	while (<IN>){
		chomp;
		next if $_ =~ /^#/;
		if ($_ =~ /^Hugo_Symbol/){
			$add_HGVSp_Short_out .= "$_\tHGVSc\tHGVSp_Short\tHGVSc\n";
		}else{
			my @line = split /\t+/,$_;
			my $k1 = "$line[3]\t$line[4]\t$line[7]";
#			print "$k1\n";
			if (exists $mut{$k1}){
				$add_HGVSp_Short_out .= "$_\t$mut{$k1}\n";
			}else{
				$add_HGVSp_Short_out .= "$_\t.\t.\t.\n";
			}
		}
	}
	close IN;
	return $add_HGVSp_Short_out;
}

sub cat_maf {
        my ($samp,$mafd) = @_;
        my %exclude = ('LOC'=>1,'ENS'=>1,'FAM'=>1,'GOL'=>1,'PRA'=>1,'NBP'=>1,'POT'=>1,'DEF'=>1,'OR2'=>1,'MUC'=>1,'KRT'=>1,'WAS'=>1,'ANK'=>1,'TRI'=>1,'OR1'=>1,'FRG'=>1,'KIAA'=>1,'CTD-'=>1,'APO'=>1);
        my $i = 0;
        my $loss = 0;
        my $maf_combine;
        foreach my $sample (@$samp){
                if (-e "$mafd/$sample.somatic.filter.maf.txt"){
                        open (IN,"$mafd/$sample.somatic.filter.maf.txt");
                        if($i == 0){
                                while (<IN>){
                                        chomp;
                                        my @tmp = split /\t/,$_;
                                        my $int = 0;
                                        foreach my $exclude(keys %exclude){
                                                if ($tmp[0] =~ /$exclude/){
                                                        $int = 1;
                                                }
                                        }
                                        next if ($int == 1);
                                        next if ($_ =~ /Splice_Region/);
                                #       next if ($_ =~ /Low_complexity/);
                                        if($tmp[4]){
                                                next if ($tmp[4] =~ /X/ || $tmp[4] =~ /Y/ || $tmp[4] =~ /chrX/ || $tmp[4] =~ /chrY/);
                                        }
                                        $maf_combine .= "$_\n";
                                }
                                $i++;
                        }else{
                                while (<IN>){
                                        chomp;
                                        my @tmp = split /\t/,$_;
                                        my $int = 0;
                                        foreach my $exclude(keys %exclude){
                                                if ($tmp[0] =~ /$exclude/){
                                                        $int = 1;
                                                }
                                        }
                                        next if ($int == 1);
                                        next if ($_ =~ /^#/);
                                        next if ($_ =~ /^Hugo_Symbol/);
                                        next if ($_ =~ /Splice_Region/);
                                #       next if ($_ =~ /Low_complexity/);
                                        if($tmp[4]){
                                                next if ($tmp[4] =~ /X/ || $tmp[4] =~ /Y/ || $tmp[4] =~ /chrX/ || $tmp[4] =~ /chrY/);
                                        }
                                        $maf_combine .= "$_\n";
                                }
                                $i++;
                        }
                        close IN;
                }else {
                       $loss++;
                        print "Warnning: can not find $mafd/$sample.somatic.filter.maf.txt,please be care!!!\n";
                }
        }
        return $maf_combine;
}

sub GetBasecounts {
	my ($samp,$bamd,$fasta,$people,$combine_maf,$od) = @_;
	my $p;
	foreach my $sample (@$samp){
		$p .= "--bam $sample:$bamd/$sample.sorted.rmdup.realigned.recal.bam ";
	}
	my $getbasecounts = "/GPFS01/softwares/GetBaseCountsMultiSample/GetBaseCountsMultiSample --omaf --fasta $fasta $p --maf $combine_maf --thread 8 --output $od/$people.GetBaseCountsMultiSample.omaf.txt --maq 10";
	return $getbasecounts;
}

sub add_chr {
	my ($maffile) = @_;
	open (IN,"$maffile");
	my $maf_add_chr;
	while (<IN>){
		chomp;
		next if $_ =~ /^#/;
		if($_ =~ /^Hugo_Symbol/){
			$maf_add_chr .= "$_\n";
		}else{
			my @all = split /\t/,$_;
			if($all[4] ne "MT"){
				$all[4] = "chr$all[4]";
			}else{
				$all[4] = "chrM";
			}
			my $tmp = join "\t",@all;
			my $tmp1 = "$tmp\n";
			$maf_add_chr .= $tmp1;
		}
	}
	close IN;
	return $maf_add_chr;
}

sub absolute {
	my ($sample,$segd,$mafd,$outd,$outd1) = @_;
        if (-e "$mafd/$sample.somatic.filter.maf.txt"){
		
		system("sed \'s/Start_Position/Start_position/g\' $mafd/$sample.somatic.filter.maf.txt >$outd/$sample.somatic.filter.maf.txt");
        }else {
		next;
        }
	if (-e "$segd/$sample.seg"){
		open (IN,"$segd/$sample.seg");
		open (OUT,">$outd/$sample.snp.seg");
		print OUT "Chromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n";
		while (<IN>){
			chomp;
			next if $_ =~ /^chrom/;
			my @all = split /\,/,$_;
			print OUT "$all[0]\t$all[9]\t$all[10]\t$all[2]\t$all[4]\n";
		}
		close OUT;
		close IN;
	}else{
		next;
	}
	my $absolute_sh = "Rscript $Bin/absolute.R $sample.somatic.filter.maf.txt $sample.snp.seg $outd $outd $outd1";
	return $absolute_sh;
}

sub Pyclone_format {
	my ($maf_file,$segd,$od) = @_;
	open (IN,"$maf_file");
	my %title;
	my $head = <IN>;
	if($head =~ /^#/){
		$head = <IN>;
	}
#	print "$head";
	die "ERROR:maf file must have title!\n" if ($head !~ /Hugo_Symbol/);
	chomp $head;
	my @title = split /\t/,$head;
	for my $i (0..$#title){
		$title{$title[$i]} = $i;
	}
	my %hash;
	while (<IN>){
		chomp;
		my @all = split /\t/,$_;
		my $k2;
		if($all[$title{HGVSp_Short}] ne "."){
			$k2 = "$all[$title{Chromosome}]--$all[$title{Start_Position}]--$all[$title{Hugo_Symbol}]--$all[$title{HGVSp_Short}]";
		}else{
			$k2 = "$all[$title{Chromosome}]--$all[$title{Start_Position}]--$all[$title{Hugo_Symbol}]--$all[$title{HGVSc}]";
		}
		my $k3 = "$all[$title{t_ref_count}]--$all[$title{t_alt_count}]--$all[$title{t_total_count}]";
		my $k1 = $all[$title{Tumor_Sample_Barcode}];
		$hash{$k1}{$k2}=$k3;
	}
	close IN;
	foreach my $k1 (keys %hash){
		foreach my $k2 (keys %{$hash{$k1}}){
			my @k2 = split /--/,$k2;
			my $format = "$hash{$k1}{$k2}--2--0--2";
			if (-e "$segd/$k1.seg"){
				open (IN,"$segd/$k1.seg");
				while (<IN>){
					chomp $_;
					next if $_ =~ /^chrom/;
					my @line = split /\,/,$_;
					$line[0] = "X" if $line[0] eq "23";
					if (($k2[0] eq $line[0]) || ($k2[0] eq "chr$line[0]")){
						if ($k2[1] >= $line[9] && $k2[1] <= $line[10]){
							my $normal_cn = 2;
							my $total_cn = $line[-2];
                                                        my $minor_cn;
                                                        if($line[-1] ne "NA"){
                                                                $minor_cn = $line[-1];
                                                        }else{
								if($total_cn == 1 || $total_cn == 0){
                                                          	     	$minor_cn = 0;
								}else{
									$minor_cn = 1;
								}
                                                        }
							my $major_cn = $total_cn - $minor_cn;
							if ($total_cn == 0){
	                                                        $major_cn = 1;
							}
                                                        $format = "$hash{$k1}{$k2}--$normal_cn--$minor_cn--$major_cn";
						}
					}else{
						next;
					}
					#last;
				}
				$hash{$k1}{$k2} = $format;
				close IN;
			}else{
				next;
			}
		}
	}
	foreach my $k1(keys %hash){
		open (FMT,">$od/$k1.PyClone.input");
		print FMT "mutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn\ttotal_depth\n";
		foreach my $k2(keys %{$hash{$k1}}){
			my @out = split /--/,$hash{$k1}{$k2};
			if($out[0] == 0 && $out[1] == 0){$out[0] = 1;$out[1] = 0;}
			print FMT "$k2\t$out[0]\t$out[1]\t$out[3]\t$out[4]\t$out[5]\t$out[2]\n";
		}
		close FMT;
	}
}

sub pid_stat{
        my $qsub_cmd = shift;
        my @cmd_id;
        if (-e "lsf_$qsub_cmd"){
                unlink glob "lsf_$qsub_cmd/* lsf_$qsub_cmd/.*";
                rmdir "lsf_$qsub_cmd" or die "cannot rmdir lsf_$qsub_cmd";
        }
        system ("/GPFS01/softwares/scripts/universal/bsubjobs_py.py -c $cpu -o $nodes $qsub_cmd");
        open CMD,'<',"./lsf_$qsub_cmd/sucesslist.csv" or die "Can't open the file lsf_$qsub_cmd/sucesslist.csv";
        while (<CMD>){
                if (/^(\d+),/){
                push @cmd_id, $1;
                }
        }
        my $a = 0 ;
        while (1){
               my $tmp1 = `bjobs`;
                my @all_id = $tmp1 =~ /^(\d+)/mg;
                my $i = 0;
                foreach my $id (@cmd_id) {
                        if (@all_id ~~ /$id$/){$i++;}
                 }
                if ($a != $i){
                        &log("$i job(s) running!");
                        $a = $i;
                }
                last  if ($i == 0);
                sleep 5;
        }
}
sub log(){
        my $message=shift;
        my $dt=DateTime->now()->set_time_zone('Asia/Chongqing');
        $message.="\t".$dt->ymd()." ".$dt->hms() ;
        print "$message\n";
}

