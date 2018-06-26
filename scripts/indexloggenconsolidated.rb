load ENV["AGILE"] + "/scripts/conf.rb"

basepath = ARGV[0]

archive = ARGV[1]

fits = Fits.new

dir = basepath + "/DATA_" + archive.to_s + "/"

fout = File.new(dir + "INDEX/LOG.log.index2", "w")

Dir[dir + "/LOG/*.gz"].sort.each do | file |
        puts file
        fits.readFitsHeader(file);
        fout.write(file + " " + fits.tt_start + " " + fits.tt_end + " LOG\n");
end
fout.close()
system("mv " + dir + "INDEX/LOG.log.index2 " + dir + "INDEX/LOG.log.index");

dir = basepath + "/FM3.119_" + archive.to_s + "/"

fout = File.new(dir + "INDEX/EVT.index2", "w")

Dir[dir + "/EVT/*.gz"].sort.each do | file |
        puts file
        if File.size(file).to_i != 0
                fits.readFitsHeader(file);
                fout.write(file + " " + fits.tt_start + " " + fits.tt_end + " EVT\n");
        end
end
fout.close()
system("mv " + dir + "INDEX/EVT.index2 " + dir + "INDEX/EVT.index")

