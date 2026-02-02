import random
import collections
import os
import zipfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def generate_dna(length=200):
	return ''.join(random.choices('ATGC', k=length))


def nucleotide_percentages(seq):
	n = len(seq)
	counts = collections.Counter(seq)
	return {nt: (counts.get(nt, 0) / n) * 100 for nt in ['A', 'T', 'G', 'C']}


def dinucleotide_percentage(seq, dinuc='CG'):
	total = max(len(seq) - 1, 1)
	count = sum(1 for i in range(len(seq) - 1) if seq[i:i+2] == dinuc)
	return (count / total) * 100, count, total


def save_results(seq, nuc_perc, cg_perc, cg_count, cg_total, out_txt='results.txt'):
	with open(out_txt, 'w') as f:
		f.write(f"Sequence ({len(seq)} bp):\n{seq}\n\n")
		f.write('Nucleotide percentages (%%):\n')
		for nt in ['A', 'T', 'G', 'C']:
			f.write(f"{nt}: {nuc_perc[nt]:.2f}\n")
		f.write('\n')
		f.write(f"Dinucleotide CG: {cg_count}/{cg_total} = {cg_perc:.4f}%%\n")


def plot_cg_percentage(cg_perc, out_img='Screenshot.jpg'):
	fig, ax = plt.subplots(figsize=(4, 6))
	ax.bar(['CG'], [cg_perc], color='#2ca02c')
	ax.set_ylim(0, 100)
	ax.set_ylabel('Percentage of CG')
	ax.set_title('CG dinucleotide percentage')
	for i, v in enumerate([cg_perc]):
		ax.text(i, v + 1, f"{v:.2f}%", ha='center')
	fig.tight_layout()
	fig.savefig(out_img, dpi=150)
	plt.close(fig)


def make_zips(workdir='.'):
	# L1.zip contains the project source and results
	l1_zip_path = os.path.join(workdir, 'L1.zip')
	with zipfile.ZipFile(l1_zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as z:
		for fname in ('L1.py', 'results.txt'):
			if os.path.exists(os.path.join(workdir, fname)):
				z.write(os.path.join(workdir, fname), arcname=fname)

	# Project_L1.zip contains ReadMe.txt, Screenshot.jpg and L1.zip
	proj_zip = os.path.join(workdir, 'Project_L1.zip')
	with zipfile.ZipFile(proj_zip, 'w', compression=zipfile.ZIP_DEFLATED) as z:
		for fname in ('ReadMe.txt', 'Screenshot.jpg', 'L1.zip'):
			if os.path.exists(os.path.join(workdir, fname)):
				z.write(os.path.join(workdir, fname), arcname=fname)


def main():
	out_dir = os.path.dirname(__file__) or '.'
	seq = generate_dna(200)
	nuc_perc = nucleotide_percentages(seq)
	cg_perc, cg_count, cg_total = dinucleotide_percentage(seq, 'CG')
	save_results(seq, nuc_perc, cg_perc, cg_count, cg_total, out_txt=os.path.join(out_dir, 'results.txt'))
	plot_cg_percentage(cg_perc, out_img=os.path.join(out_dir, 'Screenshot.jpg'))
	make_zips(out_dir)
	print('Done. Generated sequence, results.txt, Screenshot.jpg, L1.zip and Project_L1.zip')


if __name__ == '__main__':
	main()

