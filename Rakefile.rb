#!/usr/bin/ruby

require 'rake/clean'

# Set the bash environemnt variable PFUNIT
#
# export PFUNIT=/cygdrive/c/Users/John/run/pFUnit
#

LIB = "daepak-2.0.0"
STD = "f95"

MOD_DIR = File.join ".", "include"
LIB_DIR = File.join ".", "lib"
SRC_DIR = File.join ".", "src"
EX_DIR  = File.join ".", "examples"

SRC = FileList["#{SRC_DIR}/*.{f,F}*"]
OBJ = SRC.ext('o')
EX  = FileList["#{EX_DIR}/**/*.{f,F}*"]

LDFLAGS = ""
OPTIONS = "-Wextra -Wall -pedantic"

CLEAN.include("#{SRC_DIR}/*.o")
CLOBBER.include("#{LIB_DIR}/*.*",
                "#{MOD_DIR}/*.*",
                "#{EX_DIR}/**/*.o", "#{EX_DIR}/**/*.mod", "#{EX_DIR}/**/*.x")

def timestamp
  time = Time.new
  return time.strftime('%Y%m%d')
end

desc "Make the directory for modules"
directory File.basename(MOD_DIR)

desc "Make the directory for libraries"
directory File.basename(LIB_DIR)

desc "Make dependency directories in convenient task"
task :dirs do |t|
  Rake::Task[File.basename(MOD_DIR)].invoke
  Rake::Task[File.basename(LIB_DIR)].invoke
end

desc "Compile source files"
task :sources => SRC do |t|
  t.prerequisites.sort.each do |src|
    sh "gfortran -std=#{STD} -J#{MOD_DIR} #{OPTIONS} -I. -c #{src} -o #{src.ext('o')}"
  end
  sh "ar -cr #{LIB_DIR}/lib#{LIB}.a #{OBJ}"
end

desc "Run examples"
task :examples => EX do |t|
  t.prerequisites.sort.each do |src|
    begin
      dn = File.dirname(src)
      sn = FileList["#{dn}/*.f*"]
      sn.each do |s|
        puts s
        sh "gfortran -std=#{STD} #{OPTIONS} -J#{MOD_DIR} -I. -c #{s} -o #{s.ext('o')}"
        if s.include?('_prog')
          sh "gfortran -std=#{STD} -I#{dn} -I#{MOD_DIR} -L#{LIB_DIR} #{sn.ext('o')} -o #{src.ext('x')} -l#{LIB}"
          sh "#{src.ext('x')}"
        end
      end
    rescue
      puts "In rescue mode for src=[#{src}]\n\n"
      #
    end
  end
end

