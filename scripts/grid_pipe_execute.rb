#0) command to execute periodically
#1) time of sleep in seconds


begin

	a=1
	while a == 1
		cmd = ARGV[0].to_s;
		puts cmd;
		system(cmd);
		sleep (ARGV[1].to_i);
	end


end
