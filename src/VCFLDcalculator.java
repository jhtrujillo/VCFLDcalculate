import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import ngsep.variants.CalledGenomicVariant;
import ngsep.vcf.VCFFileReader;
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
	double dPrime;
	double d;

	public List<VCFRecord> loadVCF(String filename, String SNP1, String SNP2) throws IOException {
		VCFFileReader in = new VCFFileReader(filename);
		List<VCFRecord> recordsInMemory = new LinkedList<>();

		VCFRecord record1 = null;
		VCFRecord record2 = null;

		Iterator<VCFRecord> it = in.iterator();
		while (it.hasNext()) {
			VCFRecord record = it.next();
			List<CalledGenomicVariant> calls = record.getCalls();

			String snp = calls.get(0).getSequenceName() + "_" + calls.get(0).getFirst();

			if (snp.compareTo(SNP1) == 0) {
				record1 = record;
				recordsInMemory.add(record1);
			}

			if (snp.compareTo(SNP2) == 0) {
				record2 = record;
				recordsInMemory.add(record2);
				// break;
			}

		}
		return recordsInMemory;
	}

	// --------------------------------------------
	// Calculo de LD entre dos SNPs
	// --------------------------------------------
	public void CalculateLD(String filename, String SNP1, String SNP2) throws IOException {
		List<VCFRecord> recordsInMemory = loadVCF(filename, SNP1, SNP2);

		VCFRecord record1 = recordsInMemory.get(0);
		VCFRecord record2 = recordsInMemory.get(1);

		List<CalledGenomicVariant> calls1 = record1.getCalls();
		List<CalledGenomicVariant> calls2 = record2.getCalls();

		int n = calls1.size();

		int opc1 = 0;
		int opc2 = 0;
		int opc3 = 0;
		int opc4 = 0;

		// --------------------------------------------
		// Frecuencia de los alelos menores y de haplotipo AB.
		// --------------------------------------------
		for (int i = 0; i < n; i++) {
			CalledGenomicVariant call1 = calls1.get(i);
			CalledGenomicVariant call2 = calls2.get(i);

			// Dobles homocigotos.
			if (call1.isHomozygousReference() && call2.isHomozygousReference()) {
				p1 += 2;
				q1 += 2;
				frecAB += 2;
				numInd += 1;
			}

			// Homocigoto y heterocigoto
			if (call1.isHomozygousReference() && call2.isHeterozygous()) {
				p1 += 2;
				q1 += 1;
				frecAB += 1;
				numInd += 1;
			}

			// heterocigoto y homocigoto
			if (call1.isHeterozygous() && call2.isHomozygousReference()) {
				p1 += 1;
				q1 += 2;
				frecAB += 1;
				numInd += 1;
			}

		}
	}

	
	// --------------------------------------------
		// Calculo de LD entre dos SNPs
		// --------------------------------------------
		public void CalculateLDNGSEP(String filename, String SNP1, String SNP2) throws IOException {
			List<VCFRecord> recordsInMemory = loadVCF(filename, SNP1, SNP2);

			VCFRecord record1 = recordsInMemory.get(0);
			VCFRecord record2 = recordsInMemory.get(1);

			List<CalledGenomicVariant> calls1 = record1.getCalls();
			List<CalledGenomicVariant> calls2 = record2.getCalls();

			int n = calls1.size();

			int opc1 = 0;
			int opc2 = 0;
			int opc3 = 0;
			int opc4 = 0;

			// --------------------------------------------
			// Frecuencia de los alelos menores y de haplotipo AB.
			// --------------------------------------------
			for (int i = 0; i < n; i++) {
				CalledGenomicVariant call1 = calls1.get(i);
				CalledGenomicVariant call2 = calls2.get(i); 

				// Dobles homocigotos.
				if (call1.isHomozygousReference() && call2.isHomozygousReference()) {
					p1 += 2;
					q1 += 2;
					frecAB += 2;
					numInd += 1;
				}

				

			}
		}
		
		
	// --------------------------------------------
	// Calcular valores de LD NGSEP
	// --------------------------------------------
	public void getValuesLDJorge() throws IOException {
		frecAB /= (numInd * 2);
		p1 /= (numInd * 2);
		q1 /= (numInd * 2);
		p2 = 1 - p1;
		q2 = 1 - q1;
		d = frecAB - p1 * q1;

		if (p1 == 0 || q1 == 0 || p1 == 1 || q1 == 1)
			dPrime = 0;
		else if (d < 0)
			dPrime = d / Math.min(p1 * q1, (1 - p1) * (1 - q1));
		else
			dPrime = d / Math.min(p1 * (1 - q1), (1 - p1) * q1);
		r2 = d * d;
		if (p1 == 0 || q1 == 0 || p1 == 1 || q1 == 1)
			r2 = 0;
		else
			r2 /= (p1 * q1 * (1 - p1) * (1 - q1));

		// System.out.println("d\t" +"dprime\t"+ "R2");
		System.out.println(d  + "\t" + dPrime + "\t" + r2);
		
	}

	// --------------------------------------------
	// Calcular valores de LD seg�n On Measures of Gametic Disequilibrium R. C.
	// Lewontin.
	// --------------------------------------------
	public void getValuesLDLewontin() throws IOException {

		frecAB /= (numInd * 2);
		p1 /= (numInd * 2);
		q1 /= (numInd * 2);
		p2 = 1 - p1;
		q2 = 1 - q1;
		d = frecAB - p1 * q1;

		//double D = ( (p1 * p2) * (q1 * q2)) - ((p1 * p2) * (q1 * q1));
		double D = frecAB - p1 * q1;
		
		double Dmax = 0;

		if (D >= 0) {
			Dmax = Math.min(p1 * q2, q1 * p2);
			//Dmax = Math.min(p1 * (1 - q1), (1 - p1) * q1);
		} else {
			Dmax = Math.min(p1 * p2, q1 * q2);
			//Dmax = Math.min(p1 * (1 - q1), (1 - p1) * q1);
		}
		
		if (p1 == 0 || q1 == 0 || p1 == 1 || q1 == 1)
			dPrime = 0;
		else
			dPrime = D / Dmax;

		if (p1 == 0 || q1 == 0 || p1 == 1 || q1 == 1)
			r2 = 0;
		else
			r2 = (D * D) / (p1 * p2 * q1 * q2);

		// System.out.println("d\t" +"dprime\t"+ "R2");
		System.out.println(D  + "\t" + dPrime + "\t" + r2);

	}

	public static void main(String[] args) throws IOException {
		/*
		VCFLDcalculator ld1 = new VCFLDcalculator();
		ld1.CalculateLDNGSEP("/home/estuvar4/git/VCFLDcalculate/vcf/mergevcf.95ids.b_fourthfiltered-original.vcf",
				"scaffold_15457_757834", "scaffold_12739_368254");
		ld1.getValuesLDJorge();
		
		VCFLDcalculator ld2 = new VCFLDcalculator();
		ld2.CalculateLD("/home/estuvar4/git/VCFLDcalculate/vcf/mergevcf.95ids.b_fourthfiltered-original.vcf",
				"scaffold_15457_757834", "scaffold_12739_368254");
		ld2.getValuesLDLewontin();
		
		 */
		
		try {
			String opcion = args[0];

			if (opcion.compareTo("ldNGSEP") == 0 || opcion.compareTo("1") == 0) {
				try {
					VCFLDcalculator ld = new VCFLDcalculator();
					ld.CalculateLD(args[1], args[2], args[3]);
					ld.getValuesLDJorge();
					ld=null;
					// ld.getValuesLDLewontin();
				} catch (Exception e) {
					System.out.println("Try: java -jar ld.jar [ldNGSEP | 1] [path_vcf] spn1_pos spn2_pos" + e);
				}
			}

			else if (opcion.compareTo("ldLewontin") == 0 | opcion.compareTo("2") == 0) {
				try {
					VCFLDcalculator ld1 = new VCFLDcalculator();
					ld1.CalculateLD(args[1], args[2], args[3]);
					ld1.getValuesLDLewontin();
					ld1=null;
				} catch (Exception e) {
					System.out.println("Try: java -jar ld.jar [ldLewontin | 2] [path_vcf] spn1_pos spn2_pos");
				}
			}

		} catch (Exception e) {
			System.out.println("Try java -jar ld.jar -help");
		}

	}

}
