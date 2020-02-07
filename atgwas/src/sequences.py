import dataParsers, warnings
import numpy as np
import snpsdata
import pdb

"""
071009:
Resources for manipulating sequence data. E.g. handling alignments, calling SNPs, 
handling other sequence structure variations, etc.

Note that snpsdata is intended to manipulate SNP data.

Bjarni Vilhjalmsson, bvilhjal@usc.edu
"""

"""
IUPAC:
Code	 Represents		Complement  
A	Adenine			T
G	Guanine		 	C
C	Cytosine		G
T	Thymine			A
Y	Pyrimidine (C or T)	R
R	Purine (A or G)		Y
W	weak (A or T)		W
S	strong (G or C)		S
K	keto (T or G)		M
M	amino (C or A)		K
D	A, G, T (not C)		H
V	A, C, G (not T)		B
H	A, C, T (not G)		D
B	C, G, T (not A)		V
X/N	any base		X/N
-	Gap			-
"""

IUPAC_base_pair_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', '-':'-',
			'Y':'R', 'R':'Y', 'W':'W', 'S':'S', 'K':'M', 'M':'K',
			'D':'H', 'H':'D', 'V':'B', 'B':'V', 'X':'X',
			'N':'N'}


#Sorts bases by information amount. (Ad hoc) Low values are assigned to higher priority.
IUPAC_information_dict = {'A':1, 'C':1, 'G':1, 'T':1, '-':3,
			'Y':2, 'R':2, 'W':2, 'S':2, 'K':2, 'M':2,
			'D':3, 'H':3, 'V':3, 'B':3, 'X':4,
			'N':4}

class Sequence:
	"""
	A class encompassing sequence data (one sequence)
	"""
	def __init__(self, seq, seq_name, start_pos=None, end_pos=None, chr=None, quality=None, filter_threshold=10, ecotype=None):
		self.seq = seq
		self.seq_name = seq_name
		self.start_pos = start_pos
		self.end_pos = end_pos
		self.chr = chr
		self.quality = quality #Quality scores for the bases..
		if self.quality and filter_threshold > 0:
			self.filter_low_quality_calls(filter_threshold=filter_threshold)
		self.ecotype = ecotype



	def reverse_sequence(self):
		"""
		Reverses sequence, and uses complementary sequence.
		"""
		new_seq = ""
		for i in range(1, len(self.seq) + 1):
			nt = self.seq[len(self.seq) - i]
			new_seq += IUPAC_base_pair_dict[nt]
		self.seq = new_seq
		(self.start_pos, self.end_pos) = (self.end_pos, self.start_pos)


	def filter_low_quality_calls(self, filter_threshold=10, replace_char="N"):
		new_seq = ""
		for i in range(len(self.seq)):
			if self.quality[i] < filter_threshold:
				new_seq += replace_char
			else:
				new_seq += self.seq[i]
		self.seq = new_seq


	def strip(self, char='N'):
		if self.start_pos:
			i = 0
			while i < len(self.seq) and self.seq[i] == char:
				i += 1
			self.start_pos = self.start_pos + i
		self.seq = self.seq.strip(char)
		if self.end_pos:
			self.end_pos = self.start_pos + len(self.seq)

	def get_seq_diff(self, seq):
		"""
		Count's differences between this and the given sequence.
		
		Assumes they are aligned!
		"""
		if len(self.seq) != len(seq.seq):
			raise("Sequences are not equal length!")
		ignore_chars = ["N", "-"]
		diffs = 0
		counts = 0
		for i in range(len(self.seq)):
			nt1 = self.seq[i]
			nt2 = seq.seq[i]
			if not (nt1 in ignore_chars or nt2 in ignore_chars):
				counts += 1
			 	if nt1 != nt2:
					diffs += 1
		diff = -1
		if counts > 0:
			diff = diffs / float(counts)
		return {"diff":diff, "counts":counts}


	def get_sequence_boundaries(self):
		"""
		Estimates the sequence boundaries.. ignoring preceding and trailing -'s
		This is useful for handling alignments.
		
		start_pos (or self.start_pos) is assumed to start where the trailing -'s end.
		 
		"""
		i = 0
		while self.seq[i] == '-':
			i += 1
		start_i = i
		i = len(self.seq) - 1
		while self.seq[i] == '-':
			i -= 1
		end_i = i
		return (start_i, end_i)

	def update_sequence_boundaries(self):
		"""
		Replaces -'s with X's at sequence boundaries.
		"""
		i = 0
		while self.seq[i] == '-':
			i += 1
		new_seq = 'X' * (i)
		start_i = i
		i = len(self.seq) - 1
		while self.seq[i] == '-':
			i -= 1
		end_i = i + 1
		new_seq += self.seq[start_i:end_i]
		new_seq += 'X' * (len(self.seq) - i - 1)
		#print len(self.seq)
		#print len(new_seq)
		self.seq = new_seq
		return (start_i, end_i)

	def update_offsets(self, positions, offsets):
		"""
		Inserting '-' where appropriate... i.e. where there are offsets.
		"""
		seq_len = positions[-1] - positions[0]
		if seq_len != self.end_pos - self.start_pos:
			print seq_len, self.end_pos - self.start_pos
			raise Exception("lengths unequal???")
		new_seq = []
		pos_i = 0
		for i, (pos, offset) in enumerate(zip(positions, offsets)):
			if offset != 0:
				new_seq.append('-')
			else:
				new_seq.append(self.seq[pos_i])
				pos_i += 1
		new_seq = "".join(new_seq)
		self.seq = new_seq



	def __len__(self):
		"""
		Returns the sequence length ignoring -'s 
		"""
		return len(self.seq) - self.seq.count('-')


	def count_called_nts(self):
		"""
		Counts the called nucleotides, ignoring X's, N's, and -'s
		"""
		c = 0
		for nt in self.seq:
			if not nt in ["X", "N", "-"]:
				c += 1
		return c

	def merge_with(self, seq_list, priority_list=None):
		"""
		Merges this sequence with a list of sequences, which are assumed to be equally long.
		"""
		seq_len = len(self.seq)
		for s in seq_list:
			if len(s.seq) != seq_len:
				raise Exception("Sequences do not have equal length!")

		ps_list = zip(priority_list, seq_list)
		ps_list.sort()
		seq_list = map(list, zip(*ps_list))[1]
		new_seq = list(self.seq)
		for s in seq_list:
			temp_seq = list(s.seq)
			for i in range(seq_len):
				if new_seq[i] != temp_seq[i] and new_seq[i] in ['N', 'X']:
					new_seq[i] = temp_seq[i]
		self.seq = "".join(new_seq)



class SequenceData:
	"""
	A wrapper class a set of sequences..
	"""
	def __init__(self, sequences, start_pos=None, end_pos=None, chr=None, aligned=False,
			alignment_type=None, ref_seq_name=None, ref_start=None, ref_direction=None,
			ref_positions=None, ref_offsets=None):
		# start_pos end_pos should be based on relative reference sequence positions 
		self.sequences = sequences
		self.seq_names = []
		for seq in self.sequences:
			self.seq_names.append(seq.seq_name)
		#Alignment positions.. (not necessarily correct, if ref. sequence isn't aligned from the start)
		self.start_pos = start_pos
		self.end_pos = end_pos
		self.chr = chr
		self.alignment = aligned
		self.alignment_type = alignment_type
		self.ref_seq_name = ref_seq_name
		if self.ref_seq_name:
			try:
				ri = self.seq_names.index(self.ref_seq_name)
				self.ref_sequence = self.sequences[ri].seq
			except Exception, err_str:
				print self.seq_names, self.ref_seq_name
				pdb.set_trace()
				print "Unable to find the reference sequence:", self.ref_seq_name
				print err_str
				#import sys
				#sys.exit()
		self.ref_start = ref_start #When does the reference sequence start..
		self.ref_direction = ref_direction
		#print "self.start_pos,self.end_pos:",self.start_pos,self.end_pos
		#print "self.ref_start,self.ref_direction:",self.ref_start,self.ref_direction
		#print "self.ref_seq_name:",self.ref_seq_name
		if self.ref_seq_name and self.ref_start:
			self.set_reference_sequence(self.ref_seq_name, self.ref_start, ref_direction)
		elif ref_positions:
			self.ref_positions = ref_positions
		if ref_offsets:
			self.ref_offsets = ref_offsets
		self.ecotypes = []
		for seq in self.sequences:
			self.ecotypes.append(seq.ecotype)
		self.sim_m = None #Sequence similarity matrix
		self.integrity_dict = None


	def get_seq_list(self):
		"""
		Returns a list of sequences by positions, i.e. l[pos][seq]
		"""
		seq_list = []
		for seq in self.sequences:
			seq_list.append(seq.seq)
		seq_list = map(list, zip(*seq_list))
		return seq_list


	def reverse_sequences(self):
		"""
		Reverses sequences, and uses complementary sequence.
		"""
		for seq in self.sequences:
			seq.reverse_sequence()
		#(self.start_pos,self.end_pos) = (self.end_pos,self.start_pos)


	def add_data_id(self, data_id):
		for seq in self.sequences:
			seq.seq_name = data_id + "_" + seq.seq_name
		self.seq_names = []
		for seq in self.sequences:
			self.seq_names.append(seq.seq_name)


	def _guess_ecotypes_(self):
		"""
		Updates or adds ecotype IDs to data.
		"""
		import phenotypeData as pd
		a_dict = pd._getAccessionToEcotypeIdDict_(self.seq_names)
		print a_dict
		self.ecotypes = []
		for i, seq in enumerate(self.sequences):
			et = seq.ecotype
			sn = seq.seq_name.lower()
			if sn in a_dict:
				et = a_dict[sn]
			else:
				print sn, "wasn't found in dictionary."
			seq.ecotype = et
			self.ecotypes.append(et)
		print len(self.ecotypes), "ecotypes were guessed."


	def save_fasta_file(self, outfile, data_id=None):
		"""
		Saves the sequences to a file
		"""
		if data_id:
			print "Saving Fasta file with data ID:", data_id

		file_str = ""
		for i in range(len(self.seq_names)):
			#print "writing accession",self.seq_names[i]
			if data_id:
				file_str += ">" + str(data_id) + "_" + self.seq_names[i]
			else:
				file_str += ">" + self.seq_names[i]

			if self.ecotypes[i]:
				file_str += "," + str(self.ecotypes[i]) + "\n"
			else:
				file_str += "\n"

			j = -120

			for j in range(0, len(self.sequences[i].seq) - 120, 120):
				file_str += self.sequences[i].seq[j:j + 120] + "\n"
			if j < len(self.sequences[i].seq):
				file_str += self.sequences[i].seq[j + 120:] + "\n"
		f = open(outfile, "w")
		f.write(file_str)
		f.close()
		print "Saved", len(self.sequences), "sequences in the fasta file:", outfile


	def add(self, sd):
		"""
		Combines the sequence data
		"""
		if self.start_pos:
			raise NotImplementedError("Not implemented for aligned data")
		self.sequences = self.sequences + sd.sequences
		self.seq_names = self.seq_names + sd.seq_names
		self.ecotypes = self.ecotypes + sd.ecotypes


	def add_sequence(self, seq):
		self.sequences.append(seq)
		self.seq_names.append(seq.seq_name)
		self.ecotypes.append(seq.ecotype)


	def generate_unique_names(self):
		names = []
		for seq in self.sequences:
			#name = seq.seq_name.upper()
			name = seq.seq_name.lower()
			names.append(name)
			c = names.count(name)
			seq.seq_name = str(c) + "_" + name
		self.seq_names = []
		for seq in self.sequences:
			self.seq_names.append(seq.seq_name)


	def strip(self, char="N"):
		for seq in self.sequences:
			seq.strip(char)

	def filter_missing_sequences(self, filter_characters=["N", "X"]):
		new_sequences = []
		new_seq_names = []
		new_ecotypes = []
		for seq in self.sequences:
			count = 0
			for c in filter_characters:
				count += seq.seq.count(c)
			if count != len(seq.seq):
				new_sequences.append(seq)
				new_seq_names.append(seq.seq_name)
				new_ecotypes.append(seq.ecotype)
		self.sequences = new_sequences
		self.seq_names = new_seq_names
		self.ecotypes = new_ecotypes

	def filter_accessions_calls(self, missing_rate_threshold=0.25, missing_char='N'):
		new_seqs = []
		new_ecotypes = []
		new_seq_names = []
		for seq in self.sequences:
			n_count = seq.seq.count(missing_char)
			if n_count / float(len(seq.seq)) < missing_rate_threshold:
				new_seq_names.append(seq.seq_name)
				new_ecotypes.append(seq.ecotype)
				new_seqs.append(seq)
		self.sequences = new_seqs
		self.seq_names = new_seq_names
		self.ecotypes = new_ecotypes

	def get_alignment_data(self, ref_seq_name=None, ref_start=None, ref_chr=None, alignment_type=None, ref_direction=None):
		ad = AlignmentData(self.sequences, chr=ref_chr, aligned=True, alignment_type=alignment_type,
				 ref_seq_name=ref_seq_name, ref_start=ref_start, ref_direction=ref_direction)
		return ad

	def filter_missing_ecotypes(self):
		"""
		Removed accessions which don't have an ecotype.
		"""
		new_seqs = []
		new_seq_names = []
		new_ecotypes = []
		for i, e in enumerate(self.ecotypes):
			if e:
				new_seqs.append(self.sequences[i])
				new_seq_names.append(self.seq_names[i])
				new_ecotypes.append(e)
		self.seq_names = new_seq_names
		self.sequences = new_seqs
		self.ecotypes = new_ecotypes


class Indel:

	def __init__(self, start_pos, end_pos=None, length=None, type=None, lengths=None, ref_seq_name=None,
			sequences=None, alignment_pos=None, ecotypes=None):
		self.start_pos = start_pos
		if end_pos:
			self.end_pos = end_pos
		else:
			self.end_pos = start_pos #It's an insertion in the reference sequence.			
		self.length = length
		self.lengths = lengths #A list of indel lenghts for each accession.
		self.type = type
		self.ref_seq_name = ref_seq_name
		self.sequences = sequences
		self.alignment_pos = alignment_pos

	def get_binary_variant(self):
		"""
		Returns a binary variant, based on whether the insertion is in the accession or not. 
		"""
		snp = ['-'] * len(self.lengths)
		for i, l in enumerate(self.lengths):
			if l == 0:
				snp[i] = 'N'
			elif l == 'NA':
				snp[i] = 'NA'
		return snp

	def is_complex(self):
		"""
		Is this a complex indel, in the sense that it might contain more than one indels?
		"""
		return self.lengths.count(0.0) == 0

	def get_quantitative_variant(self):
		"""
		Returns a quantitative variant, based on how long the insertion is.
		"""
		raise NotImplementedError





class AlignmentData(SequenceData):
	"""
	Sequence alignment
	"""

	def _update_sequence_boundaries_(self):
		"""
		Replaces -'s with X's at sequence boundaries.
		"""
		for seq in self.sequences:
			seq.update_sequence_boundaries()

	def get_snps(self, type=0, min_called_fraction=0.1, min_indel_count=1, ignore_nts=['X', 'N'], use_letters=True, marker_type=None):
		"""
		type=0 (indels and SNPs)
		type=1 (only SNPs)
		"""
		if not self.ref_positions:
			raise Exception("A reference sequence is missing.")

		self._update_sequence_boundaries_()
		print "Retrieveing SNPs and indels from alignment."
		snps = []
		snp_positions = []
		snp_alignment_positions = []
		indels = []
		seq_list = self.get_seq_list()
		#print seq_list
		seq_num = len(self.sequences)
		min_called = min_called_fraction * seq_num
		#if len(seq_list)!=len(self.ref_positions):
		#	print len(seq_list),len(self.ref_positions)
		i = 0
		in_indel = False
		while i < len(self.ref_positions):
			#print i
			#consider this loci!
			#pos = self.ref_positions[i]
			#Does an indel begin here?
			if seq_list[i].count('-') >= min_indel_count:
				missing_count = 0
				for nt in ignore_nts:
					missing_count += seq_list[i].count(nt)

				if in_indel:#Is there a new, but adjacent indel starting?
					pass


				if missing_count < seq_num - min_called and not in_indel:
					in_indel = True
					indel_size = 0
					j = i
					indel_lengths = [0.0] * len(self.sequences)
					while j < len(self.ref_positions) and seq_list[j].count('-') >= min_indel_count: #and missing_count<seq_num-min_called
						for k in range(len(self.sequences)):
							if seq_list[j][k] == '-':
								indel_lengths[k] += 1.0
						indel_size += 1
						j += 1
					for k in range(len(self.sequences)):
						if seq_list[i][k] == 'X' and indel_lengths[k] == 0.0:
							indel_lengths[k] = 'NA'

					i_set = set(indel_lengths)
					if 'NA' in i_set:
						i_set.remove('NA')
					#else:
						#print i_set
					if len(i_set) == 1:
						print seq_list[j], j, len(self.ref_positions), len(self.sequences), len(seq_list[j])
						raise Exception
					indel = Indel(self.ref_positions[i], end_pos=self.ref_positions[j], length=j - i,
						lengths=indel_lengths, ref_seq_name=self.ref_seq_name, alignment_pos=i)
					indels.append(indel)
			elif in_indel and seq_list[i].count('-') < min_indel_count:
				in_indel = False

			nts = set(seq_list[i])
			for nt in (ignore_nts + ['-']):
				if nt in nts:
					nts.remove(nt)

			if len(nts) > 1: #It's a SNP! 
				missing_count = 0
				for nt in (ignore_nts + ['-']):
					missing_count += seq_list[i].count(nt)
				if missing_count < seq_num - min_called:
					#print "len(nts):",len(nts)
					#print "missing count:",missing_count
					#(seq_num-missing_count)>min_called

					snp = []
					for nt in seq_list[i]:
						if nt in ignore_nts:
							snp.append('NA')
						else:
							snp.append(nt)
					snps.append(snp)
					snp_positions.append(self.ref_positions[i])
					snp_alignment_positions.append(i)
			i += 1

#		print snpsd.snps[0:10]
#		print snpsd.positions[0:10]
#		print len(snpsd.snps)
#		print len(indels)
#		for i in range(10,20):
#			print indels[i].alignment_pos
#			print indels[i].start_pos
#			print indels[i].end_pos
#			print indels[i].lengths

		#print snps
		if type == 0:
			#Call the indels as snps.  Using simple 0-1..
			complex_count = 0
			for indel in indels:
				if not indel.is_complex():
					complex_count += 1
					snps.append(indel.get_binary_variant())
					snp_positions.append(indel.start_pos)
			print "Found", complex_count, "indels out of", len(indels), " which weren't complex."
			l = zip(snp_positions, snps)
			l.sort()
			l = map(list, zip(*l))
			snp_positions = l[0]
			snps = l[1]
			print "Found", len(snps), "SNPs, thereof", complex_count, "were simple indels."

		elif type == 1:
			print "Found", len(snps), "SNPs."
		else:
			raise Exception("SNP-calling type not found.")

		ecotypes = map(str, self.ecotypes)
		i_snps = []
		i_positions = []
		for indel in indels:
			i_snps.append(indel.get_binary_variant())
			i_positions.append(indel.start_pos)
		i_snpsd = snpsdata.RawSnpsData(i_snps, i_positions, accessions=ecotypes)
		#print snps[20:50]
		#print snp_positions[20:50]
		if marker_type:
			marker_types = [marker_type] * len(snps)
		else:
			marker_types = None
		snpsd = snpsdata.RawSnpsData(snps, snp_positions, accessions=ecotypes,
					alignment_positions=snp_alignment_positions, marker_types=marker_types)

		return {"snpsd":snpsd, "indels":i_snpsd}


	def get_similarity_matrix(self, indices=None):
		"""
		Direct sequence similarity (ignoring missing data, and indels) 
		"""
		if not indices:
			print "Calculating similarities between all sequences."
		else:
			print "Calculating similarities between specific sequences."
		sim_m = np.zeros((len(self.sequences), len(self.sequences)))
		count_m = np.zeros((len(self.sequences), len(self.sequences)))
		for i in range(len(self.sequences)):
			for j in range(i):
				seq1 = self.sequences[i]
				seq2 = self.sequences[j]
				if indices:
					if j in indices[i]:
						sim_m[i, j] = seq1.get_seq_diff(seq2)["diff"]
						count_m[i, j] = seq1.get_seq_diff(seq2)["counts"]
						print sim_m[i, j]
					else:
						sim_m[i, j] = None
				else:
					sim_m[i, j] = seq1.get_seq_diff(seq2)["diff"]
					count_m[i, j] = seq1.get_seq_diff(seq2)["counts"]
				sim_m[j, i] = sim_m[i, j]
				count_m[j, i] = count_m[i, j]
		self.sim_m = sim_m
		self.count_m = count_m
		return sim_m


	def check_integrity(self):
		"""
		Compares different sequences for the same ecotypes.
		"""
		e_overlaps = []
		for i, e1 in enumerate(self.ecotypes):
			overlap = []
			for j, e2 in enumerate(self.ecotypes):
				if e1 == e2 and i != j:
					overlap.append(j)
			e_overlaps.append(overlap)
		print e_overlaps
		if not self.sim_m:
			self.get_similarity_matrix(e_overlaps)

		ed = {}
		for i, e1 in enumerate(self.ecotypes):
			if len(e_overlaps[i]) > 0:
				s = set()
				for j in e_overlaps[i]:
					min_i = min(i, j)
					max_i = max(i, j)
					s.add((min_i, max_i, self.sim_m[i, j], self.count_m[i, j]))
				if e1 in ed:
					ed[e1] = ed[e1].union(s)
				else:
					ed[e1] = s

		for ei in ed:
			ed[ei] = list(ed[ei])
			print ei, ed[ei]
		self.integrity_dict = ed
		return ed


	def compare_with_snps_data(self, snpsd):
		"""
		Compares the sequence with the given (raw) snpsd for common ecotypes.
		"""
		#Ecotype map
		ecotype_map = []
		for e in self.ecotypes:
			try:
				e_i = snpsd.accessions.index(str(e))
			except Exception:
				e_i = -1
			ecotype_map.append(e_i)
		#print "ecotype map:",ecotype_map

		p_i = 0
		while snpsd.positions[p_i] < self.ref_positions[0]:
			p_i += 1


		comp_errors = []
		comp_positions = []
		comp_totals = []

		for i, pos in enumerate(self.ref_positions):
			if p_i == len(snpsd.positions):
				break
			if pos == snpsd.positions[p_i]:
				c_err = 0.0
				c_tot = 0.0
				for j, e in enumerate(self.ecotypes):
					e_i = ecotype_map[j]
					if e_i >= 0:
						nt1 = self.sequences[j].seq[i]
						nt2 = snpsd.snps[p_i][e_i]
						if not (nt1 in ['-', 'N'] or nt2 in ['-', 'N']):
							if nt1 != nt2 :
								c_err += 1.0
							c_tot += 1.0
				comp_totals.append(c_tot)
				if c_tot != 0:
					comp_errors.append(c_err / c_tot)
				else:
					comp_errors.append(-1.0)

				comp_positions.append(pos)
				p_i += 1
		print "Errors:", zip(comp_errors, comp_positions, comp_totals)
		print "Mean error:", (sum(comp_errors) / len(comp_errors))


	def merge_ecotypes(self, diff_tolerance=0.01, preserve_ref_seq=True):
		"""
		Merges common ecotypes into one sequences.		
		"""
		if not self.integrity_dict:
			self.check_integrity()
		eids = set(self.ecotypes)
		merge_eids = set()
		for ei in self.integrity_dict:
			merge_eids.add(ei)

		#Process the ecotypes to be merged
		new_seqs = []
		new_ecotypes = []
		new_seq_names = []
		for ei in merge_eids:
			diffs = {}
			counts = {}
			for (i, j, diff, count) in self.integrity_dict[ei]:
				if not i in diffs:
					diffs[i] = 0
					counts[i] = 0
				if not j in diffs:
					diffs[j] = 0
					counts[j] = 0
				if count != -1:
					diffs[i] += diff * count
					diffs[j] += diff * count
					counts[i] += count
					counts[j] += count
			for i in diffs:
				diffs[i] = diffs[i] / float(counts[i])
			for i in diffs:
				if diffs[i] > diff_tolerance:
					print self.seq_names[i], "had an error rate above the limit...but ignoring for now."
			print diffs
			temp_seqs = []
			diff_list = []
			for i in diffs:
				temp_seqs.append(self.sequences[i])
				diff_list.append(diffs[i])
			ds_list = zip(diff_list, temp_seqs)
			#print ds_list
			ds_list.sort()
			temp_seqs = []
			if preserve_ref_seq and self.sequences[self.ref_i].ecotype == ei:
				i = 0
				(d, s) = ds_list[i]
				while i < len(ds_list) - 1 and s != self.sequences[self.ref_i]:
					i += 1
					(d, s) = ds_list[i]
				(d, s) = ds_list.pop(i)
				for (d, seq) in ds_list:
					temp_seqs.append(seq)
			else:
				for (d, s) in ds_list[1:]:
					temp_seqs.append(s)
				(d, s) = ds_list[0]
			print s
			s.merge_with(temp_seqs, range(len(temp_seqs)))
			new_seqs.append(s)
			new_ecotypes.append(s.ecotype)
			new_seq_names.append(s.seq_name)

		#Add the remaining ecotypes to the new sequence list.
		remaining_eids = eids.difference(merge_eids)
		for ei in remaining_eids:
			i = self.ecotypes.index(ei)
			s = self.sequences[i]
			new_seqs.append(s)
			new_ecotypes.append(s.ecotype)
			new_seq_names.append(s.seq_name)

		self.sequences = new_seqs
		self.ecotypes = new_ecotypes
		self.seq_names = new_seq_names
		print len(self.sequences)

	def extend_alignment_with(self, ad, type=1):
		"""
		Combines the sequence data, works only on non-overlapping alignments
		
		type=1: Only common accessions
		
		type=2: All accessions 
		
		Use with care!!!
		"""
		if self.chr != ad.chr:
			raise Exception("Can't combine alignments on two different chromosomes.")
		#if self.end_pos >= ad.start_pos:
		#	raise NotImplementedError("Can't extend alignment with a preceding or overlapping alignment.")
		dist = ad.start_pos - self.end_pos - 1
		if dist > 10000:
			warnings.warn("Alignments are more more than " + str(10000) + "bp apart!")

		print "First sequence end:"
		seq_len = len(self.sequences[0].seq)
		for i in range(len(self.seq_names)):
			print self.seq_names[i], ":", self.sequences[i].seq[seq_len - 5:]
		print "Second sequence start:"
		for i in range(len(ad.seq_names)):
			print ad.seq_names[i], ":", ad.sequences[i].seq[0:5]

		print "Extending alignments, which are", dist, "apart."
		if type == 1:
			seq_map = {}
			for i, seq_name in enumerate(self.seq_names):
				if seq_name in ad.seq_names:
					seq_map[i] = ad.seq_names.index(seq_name)
			print len(self.seq_names), len(self.sequences)

			seq_list = []
			seq_names = []
			ecotypes = []
			for i in seq_map:
				seq1 = self.sequences[i]
				seq2 = ad.sequences[seq_map[i]]
				new_seq = seq1.seq + str("N" * dist) + seq2.seq
				if seq1.quality and seq2.quality:
					new_seq_quality = seq1.quality + ([0] * dist) + seq2.quality
				else:
					new_seq_quality = None
				s = Sequence(new_seq, self.seq_names[i], new_seq_quality, ecotype=self.ecotypes[i])
				seq_list.append(s)
				seq_names.append(self.seq_names[i])
				ecotypes.append(self.ecotypes[i])
		elif type == 2:
			sn = set(self.seq_names)
			new_seq_names = list(sn.union(set(ad.seq_names)))
			new_ecotypes = []
			seq_map = {}
			for i, seq in enumerate(new_seq_names):
				seq_i1 = -1
				if seq in self.seq_names:
					seq_i1 = self.seq_names.index(seq)
					ecotype = self.ecotypes[seq_i1]
				seq_i2 = -1
				if seq in ad.seq_names:
					seq_i2 = ad.seq_names.index(seq)
					ecotype = ad.ecotypes[seq_i2]
				seq_map[i] = (seq_i1, seq_i2)
				new_ecotypes.append(ecotype)

			seq_list = []
			seq_names = []
			ecotypes = []
			for i in seq_map:
				(si1, si2) = seq_map[i]
				new_seq_quality = None
				if si1 != -1 and si2 != -1:  #Both sequences exist!
					print "Both sequences exist for", self.seq_names[si1], "... merging"
					seq1 = self.sequences[si1]
					seq2 = ad.sequences[si2]
					new_seq = seq1.seq + str("N" * dist) + seq2.seq
					if seq1.quality and seq2.quality:
						new_seq_quality = seq1.quality + ([0] * dist) + seq2.quality

				elif si1 != -1:
					seq1 = self.sequences[si1]
					new_seq = seq1.seq + str("-"*(dist + len(ad.sequences[0].seq)))
				elif si2 != -1:
					seq2 = ad.sequences[si2]
					new_seq = str("-"*(dist + len(self.sequences[0].seq))) + seq2.seq
				else:
					raise Exception("Both sequences seem missing???")
				ns = Sequence(new_seq, new_seq_names[i], quality=new_seq_quality)
				seq_list.append(ns)
				seq_names.append(new_seq_names[i])
				ecotypes.append(new_ecotypes[i])
		print len(seq_map), len(seq_list)
		self.seq_names = seq_names
		self.ecotypes = ecotypes
		self.sequences = seq_list
		self.end_pos = ad.end_pos


	def get_sequence_data(self):
		"""
		Convert to SequenceData object
		
		All '-' (dashes) are removed from sequences.
		"""
		for seq in self.sequences:
			seq.seq = seq.seq.replace('-', '')
		return SequenceData(self.sequences)

	def set_reference_sequence(self, ref_seq_name, ref_start=None, ref_direction=None):
		"""
		"""
		self.ref_seq_name = ref_seq_name
		self.ref_i = self.seq_names.index(self.ref_seq_name)
		self.ref_sequence = self.sequences[self.ref_i].seq
		self.ref_start = ref_start
		self.ref_direction = ref_direction
		self._update_positions_()


	def _update_positions_(self):
		"""
		Updates the position list 'self.ref_positions'
		"""
		self.ref_positions = []
		#Find the reference sequence start position
		i = 0
		while self.ref_sequence[i] == "-":
			i += 1
		start_i = i
		for i in range(start_i):
			self.ref_positions.append(self.ref_start - start_i + i)
		pos = self.ref_start - 1
		for i in range(start_i, len(self.ref_sequence)):
			if not self.ref_sequence[i] in ['-', 'X']:
				pos += 1
			self.ref_positions.append(pos)
		self.end_pos = self.ref_positions[-1]
		if self.ref_sequence[-1] == '-':
			import warnings
			warnings.warn("Reference sequence appears to end before others end.")
		unaligned_count = self.ref_sequence.count('-') + self.ref_sequence.count('X')
		if len(self.ref_sequence) - unaligned_count != self.end_pos - self.ref_start + 1:
			import pdb
			print "Something's fishy here!"
			print unaligned_count, self.ref_start - 1
			print self.ref_positions[-1]
			print self.ref_positions[0]
			print len(self.ref_sequence) - unaligned_count, self.end_pos - self.ref_start + 1

			#pdb.set_trace()
			raise Exception



	def cut_alignment(self, start_pos, stop_pos):
		"""
		"""
		start_i = 0
		while start_pos > self.ref_positions[start_i]:
			start_i += 1
		stop_i = start_i
		while stop_i < len(self.ref_positions) and stop_pos >= self.ref_positions[stop_i]:
			stop_i += 1

		for seq in self.sequences:
			seq.seq = seq.seq[start_i:stop_i]
		self.start_pos = start_pos
		ri = self.seq_names.index(self.ref_seq_name)
		self.ref_sequence = self.sequences[ri].seq
		self._update_positions_()




	def get_allele_at_ref_pos(self, pos, window=[0, 0]):
		pos_i = self.ref_positions.index(pos)
		haplotypes = []
		for s in self.sequences:
			haplotypes.append(s.seq[pos_i - window[0]:pos_i + window[1] + 1])
		return (haplotypes, pos_i)


	def write_to_file_quan_format(self, file_name, id=None):
		out_str = ""
		if id:
			out_str += "% Alignment: " + str(id) + "\n"
		out_str += "% Chromosome:" + str(self.chr) + "\n"
		out_str += "% Start pos:" + str(self.start_pos) + "\n"
		out_str += "% End pos:" + str(self.end_pos) + "\n"
		out_str += "% Alignment length:" + str(len(self.ref_positions)) + "\n"

		#PRINT ACCESSIONS INFO
		out_str += "% Number of accessions:" + str(len(self.ecotypes)) + "\n"
		out_str += "% Ecotypes: " + (",".join(map(str, self.ecotypes))) + "\n"
		out_str += "% Position, Offset, Ecotype 1, Ecotype 2, ... \n"
		failed = False
		for i in range(len(self.ref_positions)):
			out_str += str(self.ref_positions[i]) + "," + str(self.ref_offsets[i])
			for j, seq in enumerate(self.sequences):
				if i < len(seq.seq):
					out_str += ", " + seq.seq[i]
				else:
					#print "Failed at ecotype:",self.ecotypes[j]," j =",j
					#print "Seq len:",len(seq.seq)," i =",i					 
					out_str += ", "
					failed = True
			out_str += "\n"

		if failed:
			print "Skipping sequence on chromosome", self.chr, ", start position =", self.start_pos, ", because of errors."
		else:
			f = open(file_name, "a")
			f.write(out_str + "\n\n")
			f.close()

	def print_to_screen(self, feed_num=80):
		"""
		Displays the sequence..
		"""
		pass


def readFastaFile(fastaFile, ecotype_map=None, split_name_by=None):
	"""
	Loads a Fasta file format...
	"""
	sequences = []
	f = open(fastaFile, "r")
	line = (f.readline()).strip()
	while line:
		sequence = ""
		while line[0] != ">":
			line = f.readline()
		line = line.split('>')[1]
		ecotype = None
		if ecotype_map:
			if split_name_by:
				seq_name = (line.strip()).split(split_name_by)[0]
			else:
				seq_name = line.strip()

			if seq_name in ecotype_map:
				ecotype = ecotype_map[seq_name]
			else:
				print "was unable to find", seq_name
		elif len(line) > 1:
			line = line.split(",")
			seq_name = line[0].strip()
			ecotype = int(line[1].strip())
		line = f.readline()
		while line and line[0] != ">":
			sequence += line.strip().upper()
			line = f.readline()
		sequences.append(Sequence(sequence, seq_name, ecotype=ecotype))
	f.close()
	print "Found", len(sequences), "sequences."
	num_seqs = len(sequences)
	sd = SequenceData(sequences)
	return sd



def readFastaSequence(infile, start=0, end=None, name=None, ecotype=None):
	"""
	Read a single fasta sequence...
	"""
	f = open(infile, "r")
	line = f.readline()
	while line[0] != ">":
		line = f.readline()
	line = line[1:].split(",")
	if not name:
		seq_name = line[0].strip()
	if len(line) > 1 and not ecotype:
		ecotype = int(line[1].strip())
	line = (f.readline()).strip()
	count = len(line)
	line_length = count
	while count < start:
		line = (f.readline()).strip()
		count += line_length
	cut_seq = line[start - count + line_length:]
	line = (f.readline()).strip()
	count += line_length
	if end:
		while count < end:
			cut_seq += line
			line = (f.readline()).strip()
			count += line_length
		cut_seq += line[0:end - count + line_length]
	else:
		while line:
			cut_seq += line.strip()
			line = f.readline()
	s = Sequence(cut_seq, name, ecotype=ecotype)
	return s


def readFastaAlignment(infile, ref_seq_name=None, ref_start=None, ref_chr=None, alignment_type=None, ref_direction=None):
	sd = readFastaFile(infile)
	ad = sd.get_alignment_data(ref_seq_name, ref_start, ref_chr, alignment_type, ref_direction)
	return ad


def output_phylip_format(filename, sequences, seq_names):
	#sequences = map(list,zip(*sequences))
	new_names = []
	for name in seq_names:
		if len(name) > 10:
			new_names.append(name[0:10])
		elif len(name) == 10:
			new_names.append(name)
		else:
			while len(name) < 10:
				name += " "
			new_names.append(name)


	num_seqs = len(sequences)
	num_nts = len(sequences[0])
	f = open(filename, "w")
	f.write(str(num_seqs) + " " + str(num_nts) + "\n")
	for seq, name in zip(sequences, new_names):
		seq_str = name
		for nt in seq:
			seq_str += nt
		seq_str += "\n"
		f.write(seq_str)
	f.close()




"""
def getSNPsFromFastaSequences(fastaFile,ref_start_pos=0,ref_seq=0,reversed=False):
	res = readFastaFile(fastaFile)
	sequences = res["sequences"]
	seq_names = res["names"]
#	sequences = []
#	seq_names = []
#	f = open(fastaFile,"r")
#	line = (f.readline()).strip()
#	while line:
#		sequence = ""
#		while line[0]!=">":
#			line = f.readline()
#		seq_names.append(line[1:].strip())
#		line = f.readline()
#		while line and line[0]!=">":
#			sequence += line.strip()
#			line = f.readline()
#		sequences.append(sequence)
	print "Found",len(sequences),"sequences."
	num_seqs = len(sequences)
	print seq_names 
	ref_seq = seq_names.index(ref_seq)
	print "Reference sequence was found as number:",ref_seq
	if reversed: #Then reverse sequences...
		new_sequences = []
		for seq in sequences:
			new_seq = ""
			for i in range(1,len(seq)+1):
				new_seq+=IUPAC_base_pair_dict[seq[len(seq)-i]]
			new_sequences.append(new_seq)
		
	
	sequences = map(list,zip(*new_sequences)) #transposing
	start_positions = [-1]*num_seqs
	i = 0
	while start_positions.count(-1)>2: #sequences[i].count("-")>0: #len(sequences[i])/20 or sequences[i][ref_seq]=="-":
		for j in range(0,num_seqs):
			if start_positions[j]==-1 and sequences[i][j]!="-":
				 start_positions[j] = i
		i += 1	
	start_pos = i-1
	i = 1
	seq_len = len(sequences)
	end_positions = [-1]*num_seqs
	while end_positions.count(-1)>2: #sequences[seq_len-i].count("-")>0: #len(sequences[seq_len-i])/20 or sequences[seq_len-i][ref_seq]=="-":
		for j in range(0,num_seqs):
			if end_positions[j]==-1 and sequences[seq_len-i][j]!="-":
				 end_positions[j] = seq_len-i
		i += 1
	end_pos = seq_len-i+1
	#print start_positions 
	#print end_positions

	i = start_pos
	snps_indels = []
	snps = []
	si_positions = []
	positions = []
	shifts = []
	ref_pos = start_positions[ref_seq]
	shift = 0
	print "start_pos:",start_pos
	while i < end_pos+1:
		ref_nt = sequences[i][ref_seq]
		if ref_nt!="-":
			shift = 0
			if sequences[i].count(ref_nt)!=num_seqs: #It's a SNP (or possible missing data)
				snps_indels.append(sequences[i])
				si_positions.append(ref_pos)
				if sequences[i].count("-")<1:
					snps.append(sequences[i])
					positions.append(ref_pos)
				if sequences[i].count("-")<2 and (sequences[i].count(ref_nt)+sequences[i].count("-"))<num_seqs-1:
					snps.append(sequences[i])
					positions.append(ref_pos)
				shifts.append(shift)
			ref_pos += 1
		
#		else: #reference seq has an insertion (-)
#			if sequences[i].count(ref_nt)!=num_seqs: #It's a SNP (or possible missing data)
#				snps.append(sequences[i])
#				positions.append(ref_pos)
#				shifts.append(shift)
#			shift += 1
		i += 1
	#print "ref_pos:",ref_pos
	print "end_pos:",end_pos
	new_positions = []
	for pos in positions:
		new_positions.append(ref_start_pos+start_pos+pos-1)
	positions = new_positions
	print "Number of SNPs found:",len(positions)
	return (positions,snps,seq_names)



def get_snps_from_sequences(sequences,sequence_names,ref_name=None):
	print "Found",len(sequences),"sequences."
	num_seqs = len(sequences)
	print seq_names 
	ref_seq = seq_names.index(ref_seq)
	print "Reference sequence was found as number:",ref_seq
	if reversed: #Then reverse sequences...
		new_sequences = []
		for seq in sequences:
			new_seq = ""
			for i in range(1,len(seq)+1):
				new_seq+=IUPAC_base_pair_dict[seq[len(seq)-i]]
			new_sequences.append(new_seq)
		
	
	sequences = map(list,zip(*new_sequences)) #transposing
	start_positions = [-1]*num_seqs
	i = 0
	while start_positions.count(-1)>2: #sequences[i].count("-")>0: #len(sequences[i])/20 or sequences[i][ref_seq]=="-":
		for j in range(0,num_seqs):
			if start_positions[j]==-1 and sequences[i][j]!="-":
				 start_positions[j] = i
		i += 1	
	start_pos = i-1
	i = 1
	seq_len = len(sequences)
	end_positions = [-1]*num_seqs
	while end_positions.count(-1)>2: #sequences[seq_len-i].count("-")>0: #len(sequences[seq_len-i])/20 or sequences[seq_len-i][ref_seq]=="-":
		for j in range(0,num_seqs):
			if end_positions[j]==-1 and sequences[seq_len-i][j]!="-":
				 end_positions[j] = seq_len-i
		i += 1
	end_pos = seq_len-i+1
	#print start_positions 
	#print end_positions

	i = start_pos
	snps_indels = []
	snps = []
	si_positions = []
	positions = []
	shifts = []
	ref_pos = start_positions[ref_seq]
	shift = 0
	print "start_pos:",start_pos
	while i < end_pos+1:
		ref_nt = sequences[i][ref_seq]
		if ref_nt!="-":
			shift = 0
			if sequences[i].count(ref_nt)!=num_seqs: #It's a SNP (or possible missing data)
				snps_indels.append(sequences[i])
				si_positions.append(ref_pos)
				if sequences[i].count("-")<1:
					snps.append(sequences[i])
					positions.append(ref_pos)
				if sequences[i].count("-")<2 and (sequences[i].count(ref_nt)+sequences[i].count("-"))<num_seqs-1:
					snps.append(sequences[i])
					positions.append(ref_pos)
				shifts.append(shift)
			ref_pos += 1
		
#		else: #reference seq has an insertion (-)
#			if sequences[i].count(ref_nt)!=num_seqs: #It's a SNP (or possible missing data)
#				snps.append(sequences[i])
#				positions.append(ref_pos)
#				shifts.append(shift)
#			shift += 1
		i += 1
	#print "ref_pos:",ref_pos
	print "end_pos:",end_pos
	new_positions = []
	for pos in positions:
		new_positions.append(ref_start_pos+start_pos+pos-1)
	positions = new_positions
	print "Number of SNPs found:",len(positions)
	return (positions,snps,seq_names)
	

def get_overlaping_nts(alignment_file,ref_seq,ref_start_pos,reverse=False,positions=None):
	res = readFastaFile(alignment_file)
	seqs = res["sequences"]
	seq_names = res["names"]
	#print seq_names
	#print ref_seq
	ref_index = seq_names.index(ref_seq)
	seq_eids = []
	for name in seq_names:
		seq_eids.append(FLC_dict[name])
	#print seq_eids
	if reverse:
		seqs = reverse_sequences(seqs)
	seqs = map(list,zip(*seqs)) #transposing
	ref_overlap = []
	for nts in seqs:
		if nts[ref_index]!="-":  
			ref_overlap.append(nts)
	#print len(ref_overlap),len(seqs)
	overlapping_nts = ref_overlap
	if positions:
		overlapping_nts = []
		for pos in positions:
			overlapping_nts.append(ref_overlap[pos-ref_start_pos-1])
	return (overlapping_nts,seq_eids,seq_names,ref_index)

		
def get_overlapping_snps_in_region(data_file = "/Users/bjarnivilhjalmsson/Projects/Data/250k/250K_f13_012609.csv", 
				start_pos = 3170001, end_pos = 3184001, chr=5):
	snpsd = dataParsers.parseCSVData(data_file)[chr-1]
	flc_250k_snps = []
	flc_250K_positions = []
	i = 0
	pos = snpsd.positions[i]
	while pos < start_pos:
		i+=1
		pos = snpsd.positions[i]
	while pos < end_pos:
		flc_250k_snps.append(snpsd.snps[i])
		flc_250K_positions.append(snpsd.positions[i])
		i+=1
		pos = snpsd.positions[i]
	col_250k_i = snpsd.accessions.index("6909")
	print "Col sequence from 250K:",flc_250K_positions
	#print [snp[col_250k_i] for snp in flc_250k_snps]
	
	
	(aln_snps,eids,seq_names,ref_index) = get_overlaping_nts("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/FLC_muscle_072109.aln",
								"cut_3184000_3170000",start_pos,True,flc_250K_positions)
	import analyzeHaplotype as ah
	ah.plot_haplotypes(aln_snps,accessions=seq_names,haplotypeFile="/Users/bjarnivilhjalmsson/tmp/aln_haplotype",treeFile="/Users/bjarnivilhjalmsson/tmp/aln_tree.pdf")
	eids = map(str,eids)
	print "Col sequence from FLC seq:",flc_250K_positions
	#print [snp[ref_index] for snp in aln_snps]
	
	comp_array = compareSequenceData(aln_snps,flc_250k_snps)
	print comp_array
	argmin_indices = comp_array.argmin(1)
	min_values = comp_array.min(1)
	min_eids_list = []
	for i,name in enumerate(seq_names):
		min_eids = []
		for j,val in enumerate(comp_array[i,]):
			if val <= min_values[i]:
				min_eids.append(snpsd.accessions[j])
		min_eids_list.append(min_eids)
		#print name,eids[i],(str(eids[i]) in min_eids),':', snpsd.accessions[argmin_indices[i]],':',min_eids,min_values[i]
		

	accessions = list((set(eids)).intersection(set(snpsd.accessions)))
	accessions.sort()
	print accessions
	flc_acc_map = []
	for acc in accessions:
		flc_acc_map.append(eids.index(acc))
	flc_snps = []
	for snp in aln_snps:
		new_snp = []
		for i in flc_acc_map:
			new_snp.append(snp[i])
		flc_snps.append(new_snp)
	flc_data_acc_map = flc_acc_map
	
	flc_acc_map = []
	for acc in accessions:
		flc_acc_map.append(snpsd.accessions.index(acc))
	new_flc_250k_snps = []
	for snp in flc_250k_snps:
		new_snp = []
		for i in flc_acc_map:
			new_snp.append(snp[i])
		new_flc_250k_snps.append(new_snp)
	
	flc_250k_snps = new_flc_250k_snps
	
	print flc_snps
	return (flc_250k_snps,flc_snps,flc_250K_positions,accessions,flc_data_acc_map)
	
def _comp_seqs_(seq1,seq2):
	error = 0.0
	counts = 0.0
	for nt1,nt2 in zip(seq1,seq2):
		if (not (nt1 in ['-','N'] or nt2 in ['-','N'])):
			counts += 1.0
			if nt1!=nt2:
				error +=1.0
	if counts < 5.0:
		error = 1.0
	else:
		error = error/counts
	return error

def compareSequenceData(seqs1,seqs2):
	seqs1 = map(list,zip(*seqs1))
	seqs2 = map(list,zip(*seqs2))
	import numpy as np
	comp_array = np.zeros((len(seqs1),len(seqs2)))
	for i,seq1 in enumerate(seqs1):
		for j,seq2 in enumerate(seqs2):
			comp_array[i,j] = _comp_seqs_(seq1,seq2)
	return comp_array
	


def get_overlapping_accessions_map(flc_accessions,accessions):
	flc_accs = []
	for acc in flc_accessions:
		flc_accs.append(FLC_dict[acc])
	acc_map = []
	for i in range(0,len(flc_accs)):
		if flc_accs[i] in accessions:
			acc_map.append(i)
	return acc_map
"""


def get2010sequence_datas(chr, region_start, region_end, user="bvilhjal", host="papaya.usc.edu", passwd="*rri_bjarni@usc", filter_threshold=10, id=None, ref_seq_name='Col-0'):
	"""
	Retrieve 2010 sequences within the given region.
	(Not just the SNPs)
	"""
	import dbutils
	import sys
	import time
	"""
	Retrieve 2010 data from DB.  Returns a list of RawSnpsData objects. 
	"""
	rt = time.time()
	conn = dbutils.connect_to_default_lookup()
	cursor = conn.cursor()

	#Get distinct accessions and their id.
	cmd_str = "SELECT distinct eva.ecotype_id, eva.accession_id, eva.nativename \
			FROM at.alignment an, at.accession2tg_ecotypeid eva, at.sequence s \
			WHERE an.version=3 AND s.alignment=an.id AND s.accession=eva.accession_id AND an.start<%d and an.end>%d AND an.chromosome=%d \
			ORDER BY eva.nativename" % (region_end, region_start, chr)
	acc_ecotype_dict = {}
	numRows = int(cursor.execute(cmd_str))
	accessions = []
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		e_id = int(row[0])
		a_id = int(row[1])
		a_name = row[2]
		acc_ecotype_dict[a_name] = e_id
		accessions.append(str(e_id))

	print len(accessions), " accessions found."


	#Get the distinct alignment start positions.
	cmd_str = """
	select distinct an.start, an.end, an.chromosome
	from at.alignment an, at.sequence s """
	cmd_str += "where s.alignment=an.id and an.version=3 and an.start<" + str(region_end) + " and an.end>" + str(region_start) + " and an.chromosome=" + str(chr) + " "
	cmd_str += """
	order by an.start
	"""
	numRows = int(cursor.execute(cmd_str))
	alignment_positions = []
	chromosomes = set()
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		alignment_positions.append((int(row[0]), int(row[1])))
		chromosomes.add(str(int(row[2])))

	print len(alignment_positions), " alignments found on chromosome(s)", chromosomes

	print "Fetching 2010 sequences now.."

	#Now loading all the alignments..
	sequence_datas = []
	for (start_pos, end_pos) in alignment_positions:
		print "Alignment at", start_pos, "-", end_pos
		#Offset SQL statement
		cmd_str = """
		select distinct l.position, l.offset 
		from at.alignment an, at.locus l """
		cmd_str += "where l.alignment=an.id and an.version=3 and an.start=" + str(start_pos) + " and an.end=" + str(end_pos) + " and an.chromosome=" + str(chr) + " and l.offset>0 "
		cmd_str += """
		order by l.position DESC, l.offset
		"""
		numRows = int(cursor.execute(cmd_str))
		print "	", numRows, "rows retrieved."
		off_pos_list = []
		if numRows > 0:
			curr_pos = -1
			while 1:
				row = cursor.fetchone()
				if not row:
					break;
				print row, curr_pos
				if curr_pos == -1:
					curr_pos = int(row[0])
					offset_list = [int(row[1])]
				elif curr_pos != int(row[0]):
					off_pos_list.append((curr_pos, offset_list))
					curr_pos = int(row[0])
					offset_list = [int(row[1])]
				elif curr_pos == int(row[0]):
					offset_list.append(int(row[1]))
			off_pos_list.append((curr_pos, offset_list))

		print "off_pos_list:", off_pos_list
		#Process the reference positions
		ref_positions = []
		offsets = []
		o_pos = -1
		if  len(off_pos_list):
			(o_pos, off_list) = off_pos_list.pop()
		print "o_pos:", o_pos
		for pos in range(start_pos, end_pos + 1):
			ref_positions.append(pos)
			offsets.append(0)
			if pos == o_pos:
				for o in off_list:
					ref_positions.append(pos)
					offsets.append(o)
				o_pos = -1
				if  len(off_pos_list):
					(o_pos, off_list) = off_pos_list.pop()




		#Sequence SQL statement
		cmd_str = """
		select distinct eva.nativename, s.bases, s.quality 
		from at.alignment an, at.accession2tg_ecotypeid eva, at.sequence s """
		cmd_str += "where s.alignment=an.id and s.accession=eva.accession_id and an.start=" + str(start_pos) + " and an.end=" + str(end_pos) + " and an.chromosome=" + str(chr) + " "
		cmd_str += """
		order by eva.nativename
		"""
		numRows = int(cursor.execute(cmd_str))
		print "	", numRows, "rows retrieved."
		alignment_sequences = []
		seq_names = []
		if numRows > 0:
			#i = 0
			ecotypes = []
			while(1):
				row = cursor.fetchone()
				if not row:
					break;
				#i += 1
				seq_name = row[0]
				e_i = str(acc_ecotype_dict[seq_name])
				ecotypes.append(e_i)
				seq = row[1]
				seq_quality = map(float, row[2].split(","))
				if e_i in seq_names:
					print "Ecotype already in this alignment?!"
				alignment_sequences.append(Sequence(seq, seq_name, quality=seq_quality, ecotype=e_i, filter_threshold=filter_threshold))
				#print i
			#print len(alignment_sequences)
			#print len(seq),len(ref_positions),len(offsets)
			#print ecotypes
			#print ecotypes.index('6909')
			if not '6909' in ecotypes:
				print "Reference sequence is missing... filling in the blanks with TAIR sequence."
				#col_seq = get_col_sequence(chr,start_pos,end_pos)
				#col_seq.update_offsets(ref_positions,ref_offsets)
				#alignment_sequences.append(col_seq)
				continue
			sequence_datas.append(AlignmentData(alignment_sequences, start_pos, end_pos, chr, ref_positions=ref_positions, ref_offsets=offsets, ref_seq_name=ref_seq_name))


	cursor.close ()
	conn.close ()
	dif = int(time.time() - rt)
	print "It took " + str(dif / 60) + " min. and " + str(dif % 60) + " sec. to fetch data."

	return sequence_datas

def get_col_sequence(chr, start_pos=0, end_pos=None):
	data_dir = "/Users/bjarnivilhjalmsson/Projects/Data/Col_sequence_TAIR_071009/"
	file_name = data_dir + "chr" + str(chr) + ".fas"
	return readFastaSequence(file_name, start=start_pos, end=end_pos, name="col-0", ecotype='6909')



def _get_2010_data_SNPs_n_indels_(out_file=None):
	markersds = []

	for chromosome in [1, 2, 3, 4, 5]:
		sequence_datas = get2010sequence_datas(chromosome, 0, 35000000)
		sequence_data = sequence_datas[0]
		markersd = None
		try:
			r = sequence_data.get_snps()
			markersd = r['snpsd']
			markersd.merge_data(r['indels'])
		except Exception, err_str:
			print "ERROR!", err_str

			for sequence_data in sequence_datas[1:]:
				try:
					r = sequence_data.get_snps()
					snpsd = r['snpsd']
					snpsd.merge_data(r['indels'])
					markersd.extend(snpsd)
				except Exception, err_str:
					print "ERROR!", err_str
		if markersd:
			markersds.append(markersd)
	import snpsdata
	snpsdata.writeRawSnpsDatasToFile("/tmp/2010_snps_n_indels.csv", markersds)


def _test_():
	#aln_file = "/Users/bjarnivilhjalmsson/Projects/FLC_analysis/FLC_muscle_072109.aln"
	#res = readFastaFile(aln_file)
	#(positions,aln_snps,seq_names) = getSNPsFromSequences(aln_file,ref_start_pos=3170001,ref_seq="cut_3184000_3170000",reversed=True)
	#seqs = map(list,zip(*aln_snps))
	#output_phylip_format("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/muscle_aln_snps.phylip",seqs,seq_names)
	for chrom in [1, 2, 3, 4, 5]:
		alignment_datas = get2010sequence_datas(chrom, 0, 35000000, filter_threshold=0)
		for ad in alignment_datas:
			ad.write_to_file_quan_format("/tmp/quan_2010_chr%d.csv" % chrom)



if __name__ == "__main__":
	#_run_()
	_test_()
	#_get_2010_data_SNPs_n_indels_()
	#compare_col_seqs()
	#plot_flc_haplotype()
	#print get_250k_snps_in_region()
	#print get_nts([3175000],"/Users/bjarnivilhjalmsson/Projects/FLC_analysis/FLC_muscle_072109.aln",6909,3170000,True)
	#compareRefSeqAnd250K()
	#compareFLCand250kSNPs()
	#getSNPsFromSequences("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/FLC_sequences050209.txt")
	#getSNPsFromSequences("/Users/bjarnivilhjalmsson/Projects/FLC_analysis/FLC_muscle_072109.aln")
	#cutFastaFile("/Users/bjarnivilhjalmsson/Projects/Data/Col_sequence_TAIR_071009/chr5.fas",
	#		"/Users/bjarnivilhjalmsson/Projects/Data/Col_sequence_TAIR_071009/chr5_FLC.fasta",3170000,3184000) 
