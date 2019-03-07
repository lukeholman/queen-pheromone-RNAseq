# takes a list of protein sequences and determines their orthogroup ids
#with python/2.7.8  
#takes two files, source of fasta and output
import urllib, json, time, sys
from Bio import SeqIO

urlbase = "http://www.orthodb.org/blast?"
outfile = open(sys.argv[2], "w", 0)
level = "level=7434" #Aculeata
count = 1
start = time.time()
for rec in SeqIO.parse(sys.argv[1], "fasta"):		
	while True:				
		try:
			response = urllib.urlopen(urlbase + level  + "&seq=" + str(rec.seq).replace("*",""))
			data = response.read()
			parsedData = json.loads(data)
		except:
			print "http error"
			time.sleep(2)
		else:
			break
	if parsedData['data']:
		outfile.write(rec.id + "\t" + ",".join(parsedData['data']) + "\n")
	count += 1
	if count % 100 == 0:
		print("%i %.2f" % (count, time.time() - start))
		start = time.time()		

