import numpy as np
import os
import pprint
import matplotlib.pyplot as plt

def recode(dir_path): 
	filenames = os.listdir(dir_path)
	filenames = [file for file in filenames if file[0] != '.' and file[-1] == 'z']
	filenames.append("chrX.vcf")

	for i in range(0,len(filenames) - 1):
		bash_cmd = 'vcftools --gzvcf ./Data/' + filenames[i] + ' --indv 109-EP --recode --out ./Data/' + filenames[i].split('.')[0]
		print bash_cmd
		os.system(bash_cmd)
	bash_cmd = 'vcftools --vcf ./Data/' + filenames[-1] + ' --indv 109-EP --recode --out ./Data/chrX'
	print bash_cmd
	os.system(bash_cmd)

	filenames = os.listdir(dir_path)
	filenames = [file for file in filenames if file[0] != '.' and file[5] == 'r' or file[6] == 'r']
	fw = open('./Data/indi.txt','w')
	for file in filenames:
		fr = open('./Data/' + file, 'r')
		line = fr.readline()
		while line: 
			if len(line) != 0 and line[0] != '#':
				items = line.split()
				ref_id = items[2]
				REF = items[3]
				ALT = items[4]
				gt = items[9].split(':')[0]
				print >> fw, ref_id, REF, ALT, gt
			line = fr.readline()
		fr.close()
	fw.close()

def genotype(REF, ALT, gt): 
	l, r = gt.split('|')
	l = int(l)
	r = int(r)
	if l != r: 
		z = 1
		if l == 1 and r == 0:
			geno = [ALT, REF]
			geno = "".join(geno)
		if l == 0 and r == 1:
			geno = [REF, ALT]
			geno = "".join(geno)
	if l == r:
		z = 0
		if l == 1:
			geno = [ALT, ALT]
			geno = "".join(geno)
		if l == 0:
			geno = [REF, REF]
			geno = "".join(geno)

	return z, geno

def sex(dir_path): 
	filenames = os.listdir(dir_path)
	filenames = [file for file in filenames if file[0] != '.' and file[5] == 'r' or file[6] == 'r']
	het_pcs = []
	hom_pcs = []
	for file in filenames: 
		hom = 0
		het = 0
		f = open('./Data/' + file,'r')
		line = f.readline()
		while line:
			if len(line) != 0 and line[0] != '#':
				items = line.split()
				REF = items[3]
				ALT = items[4]
				gt = items[9].split(':')[0]
				z, geno = genotype(REF, ALT, gt)
				if z == 1:
					het += 1
				if z == 0:
					hom += 1
			line = f.readline()
		het_pc = float(het)/(het + hom) * 100
		het_pcs.append(het_pc)
	plt.scatter(range(0,len(het_pcs)), het_pcs)
	plt.title('Percentages of Heterozygosity')
	plt.xlabel('Chromosome')
	plt.ylabel('Percent')
	plt.savefig('Het')

def create_categories(): 
	categories = []
	categories.append(['Blond/Red', 'Brown/Black']) # Hair color 0:0-1
	categories.append(['Straight', 'Curly']) # Hair style 1:2-3
	categories.append(['Blue/Green', 'Brown/Black']) # Eye color 2:4-5
	categories.append(['Dark', 'Light']) # Skin color 3:6-7
	categories.append(['Bitter', 'Umami']) # Taste preference 4:8-9
	categories.append(['Sense of smell']) # Sense of smell 5:10
	categories.append(['Color blindness']) # Color blindness 6:11
	cat_nums = [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], 10, 11]
	return categories, cat_nums

def create_ref_ids(): 
	ref_ids = []
	# HAIR COLOR (0:0-1)
	# Blond/Red hair
	ref_ids.append([('rs12821256', ('CC', 4), ('CT', 2), ('TT', 1)), 
		('rs12203592', ('CT', 1), ('TT', 1)), 
		('rs35264875', ('TT', 1)), 
		('rs3829241', ('AA', 5), ('AG', 2)), 
		('rs12896399', ('GG', 1)), 
		('rs1805009', ('CC', 2), ('CG', 1)), 
		('rs1805007', ('CC', 1), ('CT', 1), ('TT', 5)), 
		('rs1805008', ('CT', 1), ('TT', 5)), 
		('rs1805006', ('AA', 2), ('AC', 1)), 
		('rs11547464', ('AA', 2), ('AG', 1))])
	# Brown/Black hair
	ref_ids.append([('rs16891982', 5), 
		('rs117322171', 1), ('rs12203592', 1), ('rs12913832', 1), 
		('rs26722', 5), ('rs1667394', 1), 
		('rs183671', 3), ('rs12203592', 3), ('rs598952', 3), 
		('rs12913832', 3), ('rs1426654', 3), ('rs9820421', 1), 
		('rs58680346', 1), ('rs12355139', 1), ('rs4782309', 1), 
		('rs9949121', 1), ('rs5971087', 1)])
	# HAIR STYLE (1:2-3)
	# Straight hair
	ref_ids.append([('rs11803731', ('AA', 2)),  
		('rs7349332', ('CC', 1)), 
		('rs3827760', ('CC', 2), ('CT', 1)), 
		('rs1556547', ('AA', 1.5), ('AG', 1.75), ('GG', 2))])
	# Curly hair
	ref_ids.append([('rs11803731', ('TT', 2)),  
		('rs17646946', ('AA', 1)), 
		('rs7349332', ('TT', 2.5)), 
		('rs11803731', 3), ('rs3827760', 3), ('rs17143387', 3), 
		('rs11150606', 3), 
		('rs1556547', ('AA', 3.5), ('AG', 3.25), ('GG', 3))])
	# EYE COLOR (2:4-5)
	# Blue/Green eyes
	ref_ids.append([('rs4778138', ('AG', 1), ('GG', 1)), 
		('rs4778241', ('AA', 1), ('AC', 3), ('CC', 5)), 
		('rs7495174', ('AA', 5), ('AG', 3), ('GG', 1)), 
		('rs1129038', ('AA', 5), ('AG', 3), ('GG', 1)), 
		('rs12913832', ('GG', 5)), 
		('rs1667394', ('AA', 4)), 
		('rs1800401', ('CC', 2)), 
		('rs1800407', ('GG', 2))])
	# Brown/Black eyes
	ref_ids.append([('rs12913832', ('AA', 5), ('AG', 2)), 
		('rs1667394', ('AA', 4), ('GG', 3)), 
		('rs1800401', ('CT', 2), ('TT', 2)), 
		('rs1800407', ('AA', 3), ('AG', 3), ('CC', 2)), 
		('rs1408799', ('CT', 3)), 
		('rs12896399', ('GG', 3)), 
		('rs2733832', ('CT', 3)), 
		('rs1003719', ('AG', 3), ('CC', 3), ('GG', 3)), 
		('rs1847134', ('CC', 3), ('AC', 3)), 
		('rs1393350', ('AA', 3), ('AG', 3))])
	# SKIN COLOR (3:6-7)
	# Dark skin
	ref_ids.append([('rs26722', ('CC', 1), ('CT', 3), ('TT', 5)), 
		('rs16891982', ('CC', 5), ('CG', 3)), 
		('rs1426654', ('GG', 5))])
	# Light skin
	ref_ids.append([('rs16891982', ('GG', 3)), 
		('rs1426654', ('AA', 5), ('AG', 3)), 
		('rs1800414', 1)])
	# TASTE PREFERENCE (4:8-9)
	# Bitter
	ref_ids.append([('rs713598', ('CC', 5), ('CG', 5), ('GG', 1)), 
		('rs1726866', ('CC', 5), ('CT', 5), ('TT', 1)), 
		('rs10246939', ('CC', 5), ('CT', 5), ('TT', 1))])
	# Umami
	ref_ids.append([('rs307377', ('CT', 3), ('TT', 5))])
	# SENSE OF SMELL (5:10)
	ref_ids.append([('rs5020278', ('GG', 1), ('AG', 1), ('AA', 1)), 
		('rs6591536', ('GG', 1), ('AG', 1), ('AA', 1)), 
		('rs3751196', 3), ('rs12229599', 1), ('rs17252438', 1), 
		('rs199443', 7), ('rs6052484', 1), ('rs2251885', 2), 
		('rs2732614', 7), ('rs9321099', 3), ('rs4865875', 3), 
		('rs605843', 3), ('rs2245691', 2), ('rs2732613', 2), 
		('rs4715057', 3)])
	# COLOR BLINDNESS (6:11)
	ref_ids.append([('rs104894031', 1), 
		('rs3735972', 1), 
		('rs6471482', ('AA', 1)), 
		('rs121918344', ('CC', 1), ('CT', 1), ('TT', 1)), 
		('rs724159983', 1), ('rs104894915', 1), ('rs104894914', 1), 
		('rs104894916', 1), ('rs104894913', 1)])

	return ref_ids

def create_weights(categories): 
	weights = []
	for cat in categories:
		weights.append([0] * len(cat))
	return weights

def in_ref_ids(ref_id, geno, ref_ids):
	found_geno = 0
	matches = []

	for i in range(0,len(ref_ids)):
		for j in range(0,len(ref_ids[i])):
			if ref_id == ref_ids[i][j][0]: 
				if type(ref_ids[i][j][1]) is tuple:
					for k in range(1,len(ref_ids[i][j])):
						if geno == ref_ids[i][j][k][0]:
							found_geno = 1
							points = ref_ids[i][j][k][1]
					if found_geno:
						matches.append((i, j, k, points))
				else: 
					points = ref_ids[i][j][1]
					matches.append((i, j, points))
				found_geno = 0
			
	return matches

def char_id(match, cat_nums):
	for i in range(0,len(cat_nums)):
		if type(cat_nums[i]) is list:
			for j in range(0,len(cat_nums[i])):
				if cat_nums[i][j] == match[0]:
					cat_id = i
					subcat_id = j
		else:
			if cat_nums[i] == match[0]:
				cat_id = i
				subcat_id = 0

	return cat_id, subcat_id

def histograms(counts, weights):
	w1 = []
	w2 = []

	for i in range(0,5):
		w1.append(weights[i][0])
		w2.append(weights[i][1])
			
	fig1 = plt.figure(1)
	ax = fig1.add_subplot(111)
	indices = np.arange(5)
	width = 0.35
	ax.bar(indices-width, w1, width, color='b')
	ax.bar(indices, w2, width, color='r')

	cat_names = ['Blond/Red Hair', 'Brown/Black Hair', 
		'Straight Hair', 'Curly Hair', 
		'Blue/Green Eyes', 'Brown/Black Eyes', 
		'Dark Skin', 'Light Skin', 
		'Taste: Bitter', 'Taste: Umami']
	x = np.linspace(-0.35, 4.35, 10)
	plt.xticks(x, cat_names, rotation=25)

	rects = ax.patches
	values = []
	values.extend(w1)
	values.extend(w2)
	for rect, value in zip(rects, values):
		height = rect.get_height()
		ax.text(rect.get_x() + rect.get_width()/2, height + 1, value, 
			ha='center', va='bottom')

	plt.title('Distirbution of Weights')

	c1 = []
	c2 = []
	i = 0
	while i < 8:
		c1.append(counts[i])
		i += 1
		c2.append(counts[i])
		i += 1
	
	fig2 = plt.figure(2)
	plt.xticks([])
	plt.yticks([])
	colors = ['lightskyblue', 'lightcoral']

	ax1 = fig2.add_subplot(221)
	labels = ['Blond/Red Hair', 'Brown/Black Hair']
	sizes = [c1[0], c2[0]]
	ax1.pie(sizes, labels=labels, colors=colors, shadow=True, startangle=90, 
		autopct='%1.1f%%')

	ax2 = fig2.add_subplot(222)
	labels = ['Straight Hair', 'Curly Hair']
	sizes = [c1[1], c2[1]]
	ax2.pie(sizes, labels=labels, colors=colors, shadow=True, startangle=90, 
		autopct='%1.1f%%')

	ax3 = fig2.add_subplot(223)
	labels = ['Blue/Green Eyes', 'Brown/Black Eyes']
	sizes = [c1[2], c2[2]]
	ax3.pie(sizes, labels=labels, colors=colors, shadow=True, startangle=90, 
		autopct='%1.1f%%')

	ax4 = fig2.add_subplot(224)
	labels = ['Dark Skin', 'Light Skin']
	sizes = [c1[3], c2[3]]
	ax4.pie(sizes, labels=labels, colors=colors, shadow=True, startangle=90, 
		autopct='%1.1f%%')
	plt.title('Distribution of Counts')
	# plt.show()












# dir_path = "/Users/iammichelleau/Documents/CSE 280A/Final Project/Data"
# fr = open('./Data/indi.txt','r')

# categories, cat_nums = create_categories()
# ref_ids = create_ref_ids()
# weights = create_weights(categories)
# counts = np.zeros(len(ref_ids))

# line = fr.readline()
# while line:
# 	ref_id = line.split()[0]
# 	REF = line.split()[1]
# 	ALT = line.split()[2]
# 	gt = line.split()[3]
# 	z, geno = genotype(REF, ALT, gt)

# 	matches = in_ref_ids(ref_id, geno, ref_ids)

# 	if len(matches) > 0:
# 		for match in matches: 
# 			cat_id, subcat_id = char_id(match, cat_nums)
# 			points = match[-1]
# 			weights[cat_id][subcat_id] += points
# 			counts[match[0]] += 1
# 	line = fr.readline()
# fr.close()

# pp = pprint.PrettyPrinter(width=20)
# pp.pprint(weights)
# print counts

# fw = open('weights.txt','w')
# print >> fw, weights, counts
# fw.close()
		
# counts = [1, 11, 1, 4, 2, 3, 1, 0, 2, 0, 12, 1]
# weights = [[1, 27], [1.5, 12.5], [10, 8], [5, 0], [6, 0], [24], [1]]
histograms(counts, weights)
















