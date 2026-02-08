#this is a pipeline for the examples
#learn the format of all files by running this and check the files in "examples" directory

#run XXX.py -h for usage, some parameters could be adjusted.
#make a graph
python3 make_weighted_graph_big.py -f examples/frags_1000.txt -c examples/contacts_1000.txt -o examples/test_graph.GraphML -m ./hg19_MboI.txt

#get all top fragments (usually repeats), should be removed in the following analysis
python3 make_top_frag_dict.py -f examples/frags_1000.txt -c examples/contacts_1000.txt -t examples/top_frag_dict.npy -n 400

#save the edges of weight>2 to a dict, in order to speed up.
python3 pre_calculate_weight.py examples/test_graph.GraphML

#connection, get d-LHCCs
#python3 Components_combine_final.py -m hg19_MboI.txt -g examples/test_graph.GraphML -t examples/top_frag_dict.npy -d examples/frags_1000.txt
python3 Components_combine_RW.py -m hg19_MboI.txt -g examples/test_graph.GraphML  -t examples/top_frag_dict.npy -d examples/frags_1000.txt  -x 30 -y 14 -z 5 -i 0.08
