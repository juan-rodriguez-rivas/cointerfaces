#!/usr/bin/perl -w

use strict;
no strict "refs";

my $arg;
my $in;
my $out;
my $arg_pass        = "";
my $para            = "";
my $host_no         = 0;
my @hosts           = ();
my $cd_hit_div_exe  = "./cd-hit-div";
my $cd_hit_exe      = "./cd-hit";
my $cd_hit_2d_exe   = "./cd-hit-2d";
my $clstr_merge_exe = "./clstr_merge.pl";
my $seg_no          = 64;
my $log_file        = "";
my $restart_in      = "";
my $pbs             = 0;

my $pwd = `pwd`; chop($pwd);

while ($arg=shift) {
  if    ($arg eq "-i")       { $in         = shift; }
  elsif ($arg eq "-o")       { $out        = shift; }
  elsif ($arg eq "-B")       { $para       = shift; }
  elsif ($arg eq "-S")       { $seg_no     = shift; }
  elsif ($arg eq "-pbs")     { $pbs        = shift; }
  elsif ($arg eq "-restart") { $restart_in = shift;}
  else  {            $arg_pass .= " $arg " . shift; }
}
($in and $out) || die "not input or output";

$log_file = "$out.LOGG";

my $restart_file = "$out.restart";
my $indiv        = "$in.div";
my @commands       = ();
my @command_status = ();
my $command_no     = 0;
my $cmd;
my ($i, $j, $k, $i1, $j1, $k1);

# readin a list of hosts
if ($para) {
  open(PARA, "$para") || die "can not open $para";
  while(my $ll= <PARA>){
    chop($ll); $ll =~ s/\s//g;
    next unless ($ll);
    push(@hosts, $ll); $host_no++;
  }
  close(PARA);
}
if ($pbs) {
  for ($i=0; $i<$pbs; $i++) {
    push(@hosts, "pbs_host.$i");
  }
  $host_no = $pbs;
}
die "no host" unless $host_no;

if (-e $restart_in) {
  read_restart();
}
else {
  assign_commands();
  write_restart();
}


#dbdiv run on master node?
if ($command_status[0] eq "wait") {
  $cmd = `$commands[0]`;
  $command_status[0] = "done";
  write_restart();
}

#main runing loop
my $sleep_time = 1;
while(1) {
  #refresh job status by checking output files
  #check whether all jobs are done or not
  my $finish_flag = 1;
  my $status_change = 0;
  for ($i=1; $i<$command_no; $i++) {
    next if ($command_status[$i] eq "done");
    $finish_flag = 0;
    my $tcmd = $commands[$i];
    my $output = "";
    if ($tcmd =~ / -o\s+(\S+)/) {
      $output = $1;
      if ((-s $output) or (-s "$output.clstr")) {
        $command_status[$i] = "done";
        $status_change = 1;
      }
    }
  }
  if ($status_change) {
    write_restart();
  }
  else {
    sleep($sleep_time); print ".";
  }
  last if $finish_flag;

  my $job_sent = 0;
  for ($i=1; $i<$command_no; $i++) {
    next if ($command_status[$i] eq "done");
    next if ($command_status[$i] eq "run");
    my $tcmd = $commands[$i];
    my $in1 = "";
    my $in2 = "";
    if ($tcmd =~ / -i\s+(\S+)/) {$in1 = $1;}
    if ($tcmd =~ / -i2\s+(\S+)/) {$in2 = $1;}
    my $input_flag = 0;

    if (($in1 =~ /\S/) and ($in2 =~ /\S/)) {
      $input_flag = 1 if ((-s $in1) and (-s $in2));
    }
    elsif ($in1 =~ /\S/) {
      $input_flag = 1 if (-s $in1);
    }
    else {
      die "Error at $tcmd\n";
    }
    next unless $input_flag;

    #now input files are ready, wait
    wait_stable_file($in1);
    wait_stable_file($in2) if ($in2 =~ /\S/);

    my $thost_idx = wait_for_available_host();
    my $thost     = $hosts[$thost_idx];
    my $tsh   = "$out.$$.$thost_idx.sh";
    my $tlock = "$out.$$.$thost_idx.lock";
    my $trm   = "";
       $trm   = "rm -f $in2" if ($in2 =~ /\S/);
    open(TSH, "> $tsh") || die;
    print TSH <<EOD;
date > $tlock
$tcmd
$trm
rm -f $tlock
EOD
    close(TSH);
    if ($pbs) {
      my $t = "cd-hit-para-$thost_idx";
      open(QUEUE,"| qsub -N $t -o $t.log -e $t.err");
      print QUEUE "cd $pwd; sh $tsh";
      close(QUEUE);
    }
    else {
      $cmd = `ssh -xqf $thost 'cd $pwd; sh $tsh  >/dev/null 2>&1 &'`;
      $command_status[$i] = "run";
      print "run at $thost $tcmd\n";
    }
    $sleep_time = 1;
    $job_sent = 1;
    last;
  }

  if ((not $job_sent) and ($sleep_time < 60)) {
    $sleep_time +=5;
  }

} ############ main run loop 


######## merge all .clstr file
my $out_clstr = "$out.clstr";
if (not -s $out_clstr) {

  my @reps = ();
  for ($i=0; $i<$seg_no; $i++) {
    my $master_clstr = "$indiv-$i-o.clstr";
    die "No file $master_clstr\n" unless (-s $master_clstr);

    my $this_rep = "$indiv-$i-o";
    die "No rep $this_rep\n" unless (-e $this_rep);
    push(@reps, $this_rep);

    my @slave_clstr = ();
    for ($j=$i+1; $j<$seg_no; $j++) {
      my $tclstr = "$indiv-$j.vs.$i.clstr";
      if (-s $tclstr) {push(@slave_clstr,$tclstr); }
      else {die "No file $tclstr\n";}
    }

    if (@slave_clstr) {
      my $tclstrs = join(" ", @slave_clstr);
      print  "$clstr_merge_exe $master_clstr $tclstrs >> $out_clstr\n";
      $cmd = `$clstr_merge_exe $master_clstr $tclstrs >> $out_clstr`;
    }
    else { #this is the last piece
      print  "cat $master_clstr >> $out_clstr";
      $cmd = `cat $master_clstr >> $out_clstr`;    
    }
  }

  my $out_clstr_ren = "$out.clstr.$$";
  open(TMP, $out_clstr) || die;
  open(OTMP, "> $out_clstr_ren") || die;
  my $no = 0;
  my $cno;
  my $ll;
  while($ll=<TMP>){
    if ($ll =~ /^>Cluster (\d+)/) {
      print OTMP ">Cluster $no\n"; $no++;
      $cno  = 0;
    }
    else {
      $ll =~ s/^\d+/$cno/;
      print OTMP $ll;
      $cno++;
    }
  }
  close(TMP);
  close(OTMP);
  sleep(10);
  $cmd = `mv $out_clstr_ren $out_clstr`;

  my $reps = join(" ", @reps);
  $cmd = `cat $reps > $out`;
}


sub wait_for_available_host {
  my ($i, $j, $k);
  my $sleep = 30;
  while(1) {

    for ($i=0; $i<$host_no; $i++) {
      my $thost = $hosts[$i];
      my $tlock = "$out.$$.$i.lock";
      next if (-e $tlock);
      return $i;
    }
    sleep($sleep);
    $sleep +=30;
    if ($sleep >= 300) { $sleep = 30; }
  }
}
########## END wait_for_available_host


sub wait_stable_file {
  my ($i, $j, $k);
  my $f = shift;
  return unless (-e $f);

  my $size0 = -s $f;
  while(1) {
    sleep(10);
    my $size1 = -s $f;
    if ($size0 == $size1) { last; }
    else {$size0 = $size1; }
  }
}
########## END wait_stable_file


sub write_restart {
  my ($i, $j, $k);
  open(RES, "> $restart_file") || die;

  for ($i=0; $i<$command_no; $i++) {
    print RES "$commands[$i]\n$command_status[$i]\n";
  }
  close(RES);
}
########## END write_restart

sub assign_commands {
  my ($i, $j, $k);
  my $cmd;
  my ($idb, $idbo, $jdb, $idbout, $idblog);

  $command_no = 0;
  $cmd = "$cd_hit_div_exe -i $in -o $indiv -div $seg_no";
  push(@commands,         $cmd);
  push(@command_status, "wait");
  $command_no++;

  for ($i=0; $i<$seg_no; $i++) {
    $idb    = "$indiv-$i";
    $idblog = "$indiv-$i.log";
    #compare to previous segs
    for ($j=0; $j<$i; $j++) {
      $jdb = "$indiv-$j-o";
      $idbo = "$indiv-$i.vs.$j";
      $cmd = "$cd_hit_2d_exe -i $jdb -i2 $idb -o $idbo $arg_pass >> $idblog";
      push(@commands,         $cmd);
      push(@command_status, "wait");
      $command_no++;
      $idb = $idbo;
    }
    #self comparing
    $cmd = "$cd_hit_exe -i $idb -o $indiv-$i-o $arg_pass >> $idblog";
    push(@commands,         $cmd);
    push(@command_status, "wait");
    $command_no++;
  }
}


sub read_restart {
  $command_no = 0;
  open(RRRR, "$restart_in") || die;
  my $ll;
  while ($ll = <RRRR>) {
    chop($ll);
    push(@commands, $ll);
    $ll = <RRRR>;
    chop($ll);
    push(@command_status, $ll);
    $command_no++;
  }
  close(RRRR);
}
########## END read_restart
