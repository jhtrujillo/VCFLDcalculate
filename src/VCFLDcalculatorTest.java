import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import ngsep.variants.CalledGenomicVariant;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFLDCalculator;
import ngsep.vcf.VCFRecord;

public class VCFLDcalculator {
	double p1 = 0.0;
	double p2 = 0.0;
	double q1 = 0.0;
	double q2 = 0.0;
	double freHap = 0.0;
	double Dprime = 0;
	int numInd = 0;
	int numIndmissing = 0;
	double frecAB = 0;
	double D = 0;
	double r2 = 0;
	double r2_ngsep = 0;
	double r2_lewinton = 0;
	double dPrime;
	double d;

	public VCFRecord loadSNP(String filename, String SNP1) throws IOException {
		VCFFileReader in = new VCFFileReader(filename);

		VCFRecord recordResulto = null;

		Iterator<VCFRecord> it = in.iterator();
		while (it.hasNext()) {
			VCFRecord record = it.next();
			List<CalledGenomicVariant> calls = record.getCalls();

			String snp = calls.get(0).getSequenceName() + "_" + calls.get(0).getFirst();

			// System.out.println(snp+" "+snp.compareTo(SNP1) );

			if (snp.compareTo(SNP1) == 0) {
				// System.out.println("snp1 "+snp);
				recordResulto = record;
				break;
			}

		}

		return recordResulto;
	}

	public VCFRecord[] loadVCF(String filename, String SNP1, String SNP2) throws IOException {
		VCFFileReader in = new VCFFileReader(filename);
		VCFRecord[] recordsInMemory = new VCFRecord[2];

		recordsInMemory[0] = loadSNP(filename, SNP1);
		List<CalledGenomicVariant> calls1 = recordsInMemory[0].getCalls();

		/*
		 * for (int i = 0; i < calls1.size(); i++) {
		 * System.out.println(i+" "+calls1.get(i).isHomozygousReference()); }
		 */

		recordsInMemory[1] = loadSNP(filename, SNP2);
		return recordsInMemory;
	}

	public void CalculateLDNGSEP(String filename, String SNP1, String SNP2) throws IOException {
		VCFLDCalculator ld = new VCFLDCalculator();

		VCFRecord[] recordsInMemory = loadVCF(filename, SNP1, SNP2);

		VCFRecord record1 = recordsInMemory[0];
		VCFRecord record2 = recordsInMemory[1];

		ld.calculateLDStatistics(record1, record2);

	}

	// --------------------------------------------
	// Calcular valores de LD segun On Measures of Gametic Disequilibrium R. C.
	// Lewontin.
	// --------------------------------------------
	public void getValuesLDLewontin(String filename, String SNP1, String SNP2) throws IOException {

		VCFRecord[] recordsInMemory = loadVCF(filename, SNP1, SNP2);

		VCFRecord record1 = recordsInMemory[0];
		VCFRecord record2 = recordsInMemory[1];

		List<CalledGenomicVariant> calls1 = record1.getCalls();
		List<CalledGenomicVariant> calls2 = record2.getCalls();

		int n = calls1.size();

		p1 = 0.0;
		p2 = 0.0;
		q1 = 0.0;
		q2 = 0.0;

		int count_p1 = 0;
		int count_q1 = 0;

		freHap = 0.0;
		frecAB = 0.0;
		Dprime = 0;
		numInd = 0;
		// System.out.println("Num Ind Lewtion inicio "+ numInd);

		// System.out.println("CC n "+n);

		// --------------------------------------------
		// Frecuencia de los alelos menores y de haplotipo AB.
		// --------------------------------------------
		for (int i = 0; i < n; i++) {
			CalledGenomicVariant call1 = calls1.get(i);
			CalledGenomicVariant call2 = calls2.get(i);

			// if (call1.isUndecided() || call1.isHeterozygous()) continue;
			// if (call2.isUndecided() || call2.isHeterozygous()) continue;

			if (call1.isUndecided())
				continue;
			if (call2.isUndecided())
				continue;

			numInd++;

			// Dobles homocigotos.
			if (call1.isHomozygousReference()) {
				p1 += 1;
			}

			// Dobles homocigotos.
			if (call2.isHomozygousReference()) {
				q1 += 1;
			}

			if (call1.isHomozygousReference() && call2.isHomozygousReference()) {
				frecAB += 2;
			}

			// Homocigoto y heterocigoto
			if (call1.isHomozygousReference() && call2.isHeterozygous()) {
				frecAB += 1;
			}

			// heterocigoto y homocigoto
			if (call1.isHeterozygous() && call2.isHomozygousReference()) {
				frecAB += 1;
			}

		}

		// System.out.println("CC shared "+ numInd);

		p1 = p1 / (numInd);
		q1 = q1 / (numInd);
		p2 = 1 - p1;
		q2 = 1 - q1;
		frecAB = frecAB / (numInd * 2);

		double D = frecAB - p1 * q1;
		/*
		System.out.println("CC numInd " + numInd);
		System.out.println("CC p1 " + p1);
		System.out.println("CC q1 " + q1);
		System.out.println("CC frecAB " + frecAB);
		System.out.println("CC D " + D);
	*/
		double Dmax = 0;

		if (D >= 0) {
			Dmax = Math.min(p1 * q2, q1 * p2);
		} else {
			Dmax = Math.min(p1 * p2, q1 * q2);
		}

		if (p1 == 0 || q1 == 0 || p1 == 1.0 || q1 == 1.0)
			dPrime = 0;
		else
			dPrime = D / Dmax;

		if (p1 == 0 || q1 == 0 || p1 == 1.0 || q1 == 1.0)
			r2 = 0;
		else
			r2 = (D * D) / (p1 * p2 * q1 * q2);

		// System.out.println("d\t" +"dprime\t"+ "R2");
		System.out.println(D + "\t" + dPrime + "\t" + r2 + " " + numInd);

		this.r2_lewinton = r2;

	}

	// --------------------------------------------
	// Calcular valores de LD segun On Measures of Gametic Disequilibrium R. C.
	// Lewontin. Solo contando los sitios dobles homocigotos.
	// --------------------------------------------
	public void getValuesLDLewontinOnlyHomo(String filename, String SNP1, String SNP2) throws IOException {

		VCFRecord[] recordsInMemory = loadVCF(filename, SNP1, SNP2);

		VCFRecord record1 = recordsInMemory[0];
		VCFRecord record2 = recordsInMemory[1];

		List<CalledGenomicVariant> calls1 = record1.getCalls();
		List<CalledGenomicVariant> calls2 = record2.getCalls();

		int n = calls1.size();

		p1 = 0.0;
		p2 = 0.0;
		q1 = 0.0;
		q2 = 0.0;

		int count_p1 = 0;
		int count_q1 = 0;

		freHap = 0.0;
		frecAB = 0.0;
		Dprime = 0;
		numInd = 0;
		// System.out.println("Num Ind Lewtion inicio "+ numInd);

		// System.out.println("CC n "+n);

		// --------------------------------------------
		// Frecuencia de los alelos menores y de haplotipo AB.
		// --------------------------------------------
		for (int i = 0; i < n; i++) {
			CalledGenomicVariant call1 = calls1.get(i);
			CalledGenomicVariant call2 = calls2.get(i);

			if (call1.isUndecided() || call1.isHeterozygous())
				continue;
			if (call2.isUndecided() || call2.isHeterozygous())
				continue;

			numInd++;

			// Dobles homocigotos.
			if (call1.isHomozygousReference()) {
				p1 += 1;
			}

			// Dobles homocigotos.
			if (call2.isHomozygousReference()) {
				q1 += 1;
			}

			if (call1.isHomozygousReference() && call2.isHomozygousReference()) {
				frecAB += 2;
			}

		}

		// System.out.println("CC shared "+ numInd);

		p1 = p1 / (numInd);
		q1 = q1 / (numInd);
		p2 = 1 - p1;
		q2 = 1 - q1;
		frecAB = frecAB / (numInd * 2);

		double D = frecAB - p1 * q1;
		/*
		System.out.println("CC numInd " + numInd);
		System.out.println("CC p1 " + p1);
		System.out.println("CC q1 " + q1);
		System.out.println("CC frecAB " + frecAB);
		System.out.println("CC D " + D);
		 */
		double Dmax = 0;

		if (D >= 0) {
			Dmax = Math.min(p1 * q2, q1 * p2);
		} else {
			Dmax = Math.min(p1 * p2, q1 * q2);
		}

		if (p1 == 0 || q1 == 0 || p1 == 1.0 || q1 == 1.0)
			dPrime = 0;
		else
			dPrime = D / Dmax;

		if (p1 == 0 || q1 == 0 || p1 == 1.0 || q1 == 1.0)
			r2 = 0;
		else
			r2 = (D * D) / (p1 * p2 * q1 * q2);

		// System.out.println("d\t" +"dprime\t"+ "R2");
		System.out.println(D + "\t" + dPrime + "\t" + r2 + " " + numInd);

		this.r2_lewinton = r2;

	}

	public static void main(String[] args) throws IOException {
		/*
		VCFLDcalculator ld1 = new VCFLDcalculator();

		String snp1 = "Pl01_";
		String snp2 = "Pl01_";

		// snp1 += "122138";
		// snp2 += "122025";

		snp1 += "274067";
		snp2 += "267787";

		// snp1 += "273978";
		// snp2 += "122025";

		ld1.CalculateLDNGSEP("/home/estuvar4/Desktop/tmp.vcf", snp1, snp2);
		System.out.println();
		ld1.getValuesLDLewontin("/home/estuvar4/Desktop/tmp.vcf", snp1, snp2);
		System.out.println();
		ld1.getValuesLDLewontinOnlyHomo("/home/estuvar4/Desktop/tmp.vcf", snp1, snp2);
		*/

		try {
			String opcion = args[0];

			if (opcion.compareTo("ldNGSEP") == 0 || opcion.compareTo("1") == 0) {
				try {
					VCFLDcalculator ld = new VCFLDcalculator();
					ld.CalculateLDNGSEP(args[1], args[2], args[3]);
					ld=null;
				} catch (Exception e) {
					System.out.println("Try: java -jar ld.jar [ldNGSEP | 1] [path_vcf] spn1_pos spn2_pos" + e);
				}
			}

			else if (opcion.compareTo("ldLewontin") == 0 | opcion.compareTo("2") == 0) {
				try {
					VCFLDcalculator ld = new VCFLDcalculator();
					ld.getValuesLDLewontin(args[1], args[2], args[3]);
					ld=null;
				} catch (Exception e) {
					System.out.println("Try: java -jar ld.jar [ldLewontin | 2] [path_vcf] spn1_pos spn2_pos | " +e);
				}
			}
			
			else if (opcion.compareTo("ldLewontinHomo") == 0 | opcion.compareTo("3") == 0) {
				try {
					VCFLDcalculator ld = new VCFLDcalculator();
					ld.getValuesLDLewontinOnlyHomo(args[1], args[2], args[3]);
					ld=null;
				} catch (Exception e) {
					System.out.println("Try: java -jar ld.jar [ldLewontinHomo | 3] [path_vcf] spn1_pos spn2_pos | " +e);
				}
			}
		
		} catch (Exception e) {
			System.out.println("Try java -jar ld.jar -help");
		}
		
		 
		
	}

}
