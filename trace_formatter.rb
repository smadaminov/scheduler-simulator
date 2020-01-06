count = 0
found = 0
num_active_jobs = 0

folder_dir = ARGV[0]
folder_dir = './' if ARGV.empty?

# Formatting file with start nodes

file_src = folder_dir + 'startingNodes.csv'
file_dst = 'startingNodes.dat'
filename = folder_dir + 'prednodes.csv'

if not File.file?(file_src)
    puts 'Cannot find file with start nodes (startingNodes.csv).'
    exit
end

if not File.file?(filename)
    puts 'Cannot find file with predicate nodes (prednodes.csv).'
    exit
end

puts "Parsing #{file_src}..."

# Let's first create a dictionary to store correspondence
# between ID and type

types = {}
File.readlines(filename).drop(1).each do |str|
    id, type = str.delete("\n").delete("\"").split("|")
    types[type] = id
end
puts "There are #{types.length} types"

File.open(file_dst, "w") do |file|
    File.readlines(file_src).drop(1).each do |line|
        if types.key? line.delete!("\n")
            puts types[line]
            file.puts "#{types[line]}"
            found = found + 1
        else
            #p "Found non-existing type."
            #exit
        end
        count = count + 1
    end
end

file_src = folder_dir + 'startingNodes.csv'
file_dst = 'startingNodes.dat'
filename = folder_dir + 'unitnodes.csv'

types = {}
File.readlines(filename).drop(1).each do |str|
    id, type = str.delete("\n").delete("\"").split("|")
    types[type] = id
end
puts "There are #{types.length} types"

File.open(file_dst, "a") do |file|
    File.readlines(file_src).drop(1).each do |line|
        if types.key? line.delete!("\n")
            puts types[line]
            file.puts "#{types[line]}"
            found = found + 1
        else
            #p "Found non-existing type."
            #p types[line]
            #exit
        end
        count = count + 1
    end
end

puts "Done parsing file with start nodes."

# Formatting file with processing times

file_src = folder_dir + "processingTime.csv"
file_dst = 'processingTime.dat'

if not File.file?(file_src)
    puts "Cannot find file with processing times (processingTime.csv)."
    exit
end

puts "Parsing #{file_src}..."

File.open(file_dst, "w") do |file|
    File.readlines(file_src).drop(1).each do |line|
        id, time = line.delete("\n").split("|")
        file.puts "#{id} #{time}"
        num_active_jobs = num_active_jobs + 1
    end
end

file_dst = 'activatedJob.dat'

File.open(file_dst, "w") do |file|
    File.readlines(file_src).drop(1).each do |line|
        id, time = line.delete("\n").split("|")
        file.puts "#{id}"
    end
end

puts "Done parsing file with processing time."

puts "Found #{found}/#{count} start nodes."
puts "There #{num_active_jobs} active jobs."

puts "Done."
