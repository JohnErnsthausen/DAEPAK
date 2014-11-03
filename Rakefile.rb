#!/usr/bin/ruby

require 'rake/clean'

LIB = "daepak"

ATOM = File.join ".", "atom"
SRC_DIR = File.join ".", "src"
EX_DIR  = File.join ".", "examples"
EXAMPLES = FileList["#{EX_DIR}/*/"].uniq.sort

OPTIONS = "-Wall -pedantic"

CLEAN.include("#{ATOM}","./**/*.{o,mod,x}")
CLOBBER.include("./**/*.{log,txt}", "lib#{LIB}.a")

desc "Compile source files"
task :sources do |t|
  sh "cp -r #{SRC_DIR} #{ATOM}"
  Dir.chdir(ATOM) do
    FileList["#{"./*.{f,F}*"}"].uniq.sort.each do |src|
      # Remove the prefix on the source file that orders compilation
      srcfile = File.basename(src)
      outfile = File.basename(src).ext('o')
      sh "gfortran #{OPTIONS} -c #{srcfile} -o #{outfile}"
    end
  end
  sh "ar -cr lib#{LIB}.a #{FileList["#{ATOM}/**/*.{o,mod}"]}"
end

desc "Run examples"
task :examples => EXAMPLES do |t|
  t.prerequisites.uniq.sort.each do |ex|
    # Order is important here. Get example sources. Then copy in solver sources.
    puts sources = FileList["#{ex}/*.f*"].collect{|f| File.basename(f)}
    sh "cp lib#{LIB}.a #{ex}"
    Dir.chdir(ex) do
      begin
        sh "ar x lib#{LIB}.a && rm lib#{LIB}.a"
        sources.each do |s|
          puts s
          sh "gfortran #{OPTIONS} -c #{s} -o #{s.ext('o')}"
          if s.include?('_prog')
            puts objects = FileList["*.o"].join(' ')
            sh "gfortran #{OPTIONS} #{objects} -o #{s.ext('x')}"
            sh "./#{s.ext('x')}"
          end
        end
      rescue
        puts "In rescue mode for src=[#{src}]\n\n"
        #
      end
    end
  end
end

desc "Run a single example"
task :example, :ex do |t, args|
  model = "#{EX_DIR}/#{args[:ex]}/"
  raise unless EXAMPLES.include(model)
  # Order is important here. Get example sources. Then copy in solver sources.
  puts sources = FileList["#{model}/*.f*"].collect{|f| File.basename(f)}
  sh "cp lib#{LIB}.a #{model}"
  Dir.chdir(model) do
    begin
      sh "ar x lib#{LIB}.a && rm lib#{LIB}.a"
      sources.each do |s|
        puts s
        sh "gfortran #{OPTIONS} -c #{s} -o #{s.ext('o')}"
        if s.include?('_prog')
          puts objects = FileList["*.o"].join(' ')
          sh "gfortran #{OPTIONS} #{objects} -o #{s.ext('x')}"
          sh "./#{s.ext('x')}"
        end
      end
    rescue
      puts "In rescue mode for src=[#{model}]\n\n"
      #
    end
  end
end




