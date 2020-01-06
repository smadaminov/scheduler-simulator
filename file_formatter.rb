require 'optparse'

dict = {}
count = 0
seq_count = 0
num_edges = 0
num_prednodes = 0

ARGV << '-h' if ARGV.empty?

options = {}
OptionParser.new do |opt|
    opt.on('-e', '--edges=EDGES_FILE', 'File with edges') { |o| options[:edges_file] = o }
    opt.on('-p', '--predicates=PREDICATES_FILE', 'File with predicate nodes') { |o| options[:predicates_file] = o }
    opt.on('-s', '--sequence=SEQUENCE_FILE', 'Sequence file with node priorities') { |o| options[:sequence_file] = o }
    opt.on_tail('-h', '--help', 'Show help message') do
        puts opt
        exit
    end
end.parse!

edges_file = options[:edges_file]
predicates_file = options[:predicates_file]
sequence_file = options[:sequence_file]

if not File.file?(edges_file)
    puts "Please provide correct file name for file with edges."
    exit
end

if not File.file?(predicates_file)
    puts "Please provide correct file name for file with predicates."
    exit
end

if not File.file?(sequence_file)
    puts "Please provide correct file name for file with sequence numbers."
    exit
end

puts "Parsing #{edges_file}..."

File.readlines(edges_file).drop(1).each do |line|
    src, dst = line.delete!("\n").split('|')
    if not dict.key?(src)
        dict[src] = count
        count = count + 1
    end
    if not dict.key?(dst)
        dict[dst] = count
        count = count + 1
    end
    num_edges = num_edges + 1
end

file_dst = 'edges.dat'
File.open(file_dst, 'w') do |file|
    File.readlines(edges_file).drop(1).each do |line|
        src, dst = line.delete!("\n").split('|')
        file.puts "#{src} #{dst} #{dict[src]} #{dict[dst]}"
    end
end

puts "Parsing #{sequence_file}..."

file_dst = 'sequence.dat'
File.open(file_dst, 'w') do |file|
    File.readlines(sequence_file).drop(1).each do |line|
        seq_count = seq_count + 1
        node, priority = line.delete!("\n").split('|')
        if not dict.key?(node)
            puts "ERROR: Unknown node. Exiting..."
            exit
        end
        file.puts "#{node} #{dict[node]} #{priority}"
    end
end

if not seq_count.eql? count
    puts "ERROR: Mismatch in number of sequence priorities and number of nodes. Exiting..."
    exit
end

puts "Parsing #{predicates_file}..."

file_dst = "prednodes.dat"
File.open(file_dst, 'w') do |file|
    File.readlines(predicates_file).drop(1).each do |line|
        prednode = line.partition('|').first
        file.puts prednode
        num_prednodes = num_prednodes + 1
    end
end


puts "There are #{count} unique vertices."
puts "There are #{num_edges} edges."
puts "There are #{num_prednodes} predicate nodes."

puts 'Done.'
