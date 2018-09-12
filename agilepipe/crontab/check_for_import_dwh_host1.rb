command_name = "import_dwh_host1.sh"

  job_list =  %x[squeue --format=%j]

  job_list =  job_list.split("\n")
  find = false
  job_list.each do |job_name|
    #puts line
    if(job_name.delete("\n") == command_name )
      puts "gia esistente"
      find = true
    end

  end

if(!find)
  puts "lancio comando"
  system("sbatch --partition=large --output=/tmp/import_dwh.out /opt/prod/AGILEPIPE/crontab/"+command_name+" ")
end
