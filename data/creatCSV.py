out = open("chr-22.geno.reduced.csv","w")
with open("chr-22.geno.reduced") as f:
	for line in f:
		out.write(",".join(list(line.strip()))+"\n")

out.close()