out = open("chr-20.geno.csv","w")
with open("chr-20.geno") as f:
	for line in f:
		out.write(",".join(list(line.strip()))+"\n")

out.close()