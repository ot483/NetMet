# NetMet
NetMet: A Network-Based Tool for Predicting Metabolic Capacities of Microbial Species and their Interactions <br />

This is the local version of the tool published at: https://pubmed.ncbi.nlm.nih.gov/32503277/ <br /> authors: 
Ofir Tal, Gopinath Selvaraj, Shlomit Medina, Shany Ofaim and Shiri Freilich




# Usage: 
python NetMet.py <BaseFolder path> <genomes file> <envs file>

# Example:
  
python NetMet.py ./ ./genomes.txt ./envs.txt


# env file format: 

Index|Env_Content <br />
1|4: 87 148 112 53 <br />
2|2: 87 148 <br />
3|4: 87 148 12 2 <br />



    

# genomes file format:

Ham_Rao_cleaned 6.4.1.2 3.1.1.85 2.3.1.47 2.6.1.1 1.-.-.- 1.1.1.100 1.1.1.193 1.1.1.205 1.1.1.22 1.1.1.25 1.1.1.262 1.1.1.271 .... <br />
Rick_cleaned 5.4.99.12 4.2.1.46 4.3.1.19 2.3.1.180 3.4.21.53 2.5.1.7 2.7.4.8 4.1.99.3 2.7.9.1 .... <br />
 ... 
