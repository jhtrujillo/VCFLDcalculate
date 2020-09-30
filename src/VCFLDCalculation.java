import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;
import ngsep.variants.VariantCallReport;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFLDCalculator;
import ngsep.vcf.VCFRecord;

public class VCFLDCalculation {
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

	// --------------------------------------------
	// Calcular valores de LD segun On Measures of Gametic Disequilibrium R. C.
	// Lewontin.
	// --------------------------------------------
	public void getpairseLDLewontin(VCFRecord record1, VCFRecord record2) throws IOException {

		List<CalledGenomicVariant> calls1 = record1.getCalls();
		List<CalledGenomicVariant> calls2 = record2.getCalls();

		String SNP1 = calls1.get(0).getSequenceName() + " " + calls1.get(0).getFirst();
		String SNP2 = calls2.get(0).getSequenceName() + " " + calls2.get(0).getFirst();

		int n = calls1.size();

		p1 = 0.0;
		p2 = 0.0;
		q1 = 0.0;
		q2 = 0.0;

		freHap = 0.0;
		frecAB = 0.0;
		Dprime = 0;
		numInd = 0;

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

		p1 = p1 / (numInd);
		q1 = q1 / (numInd);
		p2 = 1 - p1;
		q2 = 1 - q1;
		frecAB = frecAB / (numInd * 2);

		double D = frecAB - p1 * q1;

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

		if (dPrime > 1)
			dPrime = 1;
		
		if (r2 > 1)
			r2 = 1;

		System.out.println(SNP1 + "\t" + SNP2 + "\t" + D + "\t" + dPrime + "\t" + r2 + "\t" + numInd);

		this.r2_lewinton = r2;

	}

	// --------------------------------------------
	// Calcular valores de LD segun On Measures of Gametic Disequilibrium R. C.
	// Lewontin. Solo contando los sitios dobles homocigotos.
	// --------------------------------------------
	public void getValuesLDLewontinOnlyHomo(VCFRecord record1, VCFRecord record2) throws IOException {

		List<CalledGenomicVariant> calls1 = record1.getCalls();
		List<CalledGenomicVariant> calls2 = record2.getCalls();

		String SNP1 = calls1.get(0).getSequenceName() + " " + calls1.get(0).getFirst();
		String SNP2 = calls2.get(0).getSequenceName() + " " + calls2.get(0).getFirst();

		int n = calls1.size();

		p1 = 0.0;
		p2 = 0.0;
		q1 = 0.0;
		q2 = 0.0;

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

			if (call1.isUndecided() || call1.isHeterozygous())	continue;
			if (call2.isUndecided() || call2.isHeterozygous()) continue;

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
		 * System.out.println("CC numInd " + numInd); System.out.println("CC p1 " + p1);
		 * System.out.println("CC q1 " + q1); System.out.println("CC frecAB " + frecAB);
		 * System.out.println("CC D " + D);
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
		
		
		if (dPrime > 1)
			dPrime = 1;
		
		if (r2 > 1)
			r2 = 1;

		System.out.println(SNP1 + "\t" + SNP2 + "\t" + D + "\t" + dPrime + "\t" + r2 + "\t" + numInd);

		this.r2_lewinton = r2;

	}

	
	public void CalculateLDNGSEP(VCFRecord record1, VCFRecord record2) throws IOException {
		VCFLDCalculator ld = new VCFLDCalculator();
		ld.calculateLDStatistics(record1, record2);

	}

	public void VCFldcalulation(String filename, String Chr, String pos, String opcion) throws IOException {

		String SNP1 = Chr + "_" + pos;

		VCFFileReader in = new VCFFileReader(filename);
		VCFRecord SNPrecord1 = null;
		Iterator<VCFRecord> it = in.iterator();

		// Saco el VCFRecord del SNP de interes.
		while (it.hasNext()) {
			VCFRecord record = it.next();
			List<CalledGenomicVariant> calls = record.getCalls();

			String snp = calls.get(0).getSequenceName() + "_" + calls.get(0).getFirst();

			if (snp.compareTo(SNP1) == 0) {
				SNPrecord1 = record;
			}

		}

			
		// Comparo el VCFRecord seleccionado con el resto de SNPs.
		VCFFileReader in2 = new VCFFileReader(filename);
		Iterator<VCFRecord> it2 = in2.iterator();
		System.out.println("SNP1" + "\t" + "SNP2" + "\t" + "D" + "\t" + "dPrime" + "\t" + "r2" + "\t" + "numInd");

		while (it2.hasNext()) {
			VCFRecord SNPrecord2 = it2.next();

			if (opcion.compareTo("Lewontin") == 0) {
				getpairseLDLewontin(SNPrecord1, SNPrecord2);
			} else if (opcion.compareTo("LewontinOnlyHomo") == 0) {
				getValuesLDLewontinOnlyHomo(SNPrecord1, SNPrecord2);
			} 
			else if (opcion.compareTo("NGSEP") == 0) {
				CalculateLDNGSEP(SNPrecord1, SNPrecord2);
			}

		}

	}

	public void VCFldcalulationAll(String filename, String opcion) throws IOException {

		VCFFileReader in = new VCFFileReader(filename);
		Iterator<VCFRecord> it = in.iterator();

		System.out.println("SNP1" + "\t" + "SNP2" + "\t" + "D" + "\t" + "dPrime" + "\t" + "r2" + "\t" + "numInd");

		// Saco el VCFRecord del SNP de interes.
		while (it.hasNext()) {
			VCFRecord SNPrecord1 = it.next();

			// Comparo el VCFRecord seleccionado con el resto de SNPs.
			VCFFileReader in2 = new VCFFileReader(filename);
			Iterator<VCFRecord> it2 = in2.iterator();

			while (it2.hasNext()) {
				VCFRecord SNPrecord2 = it2.next();

				if (opcion.compareTo("Lewontin") == 0) {
					getpairseLDLewontin(SNPrecord1, SNPrecord2);
				} else if (opcion.compareTo("LewontinOnlyHomo") == 0) {
					getValuesLDLewontinOnlyHomo(SNPrecord1, SNPrecord2);
				} else if (opcion.compareTo("NGSEP") == 0) {
					CalculateLDNGSEP(SNPrecord1, SNPrecord2);
				}
			}
		}
	}

	public static void main(String[] args) throws IOException {

		//VCFLDCalculation ld = new VCFLDCalculation();

		//String Chr = "Pl01";
		//String snp = "121870";
		
		//ld.VCFldcalulation("/home/estuvar4/Desktop/tmp.vcf", Chr, snp, "LewontinOnlyHomo");

		// ld.VCFldcalulation("/home/estuvar4/Desktop/tmp.vcf", Chr, snp, "Lewontin");

		// ld.VCFldcalulation("/home/estuvar4/Desktop/tmp.vcf", Chr, snp,
		// "LewontinOnlyHomo");

		// ld.VCFldcalulation("/home/estuvar4/Desktop/tmp.vcf", Chr, snp, "NGSEP");
		
		//ld.VCFldcalulationAll("/home/estuvar4/Desktop/tmp.vcf", "Lewontin");
		//ld.VCFldcalulationAll("/home/estuvar4/Desktop/tmp.vcf", "LewontinOnlyHomo");
		//ld.VCFldcalulationAll("/home/estuvar4/Desktop/tmp.vcf", "NGSEP");
		
		
		
		try {
			String opcion = args[0];
			
			if (opcion.compareTo("-h") == 0 || opcion.compareTo("-help") == 0) {
				try {
					System.out.println("Try: java -jar ld.jar [VCFldcalculator | 1] [path_vcf] Method[Lewontin|LewontinOnlyHomo|NGSEP]");
					System.out.println("Try: java -jar ld.jar [VCFldcalculatorPairse | 2] [path_vcf] Chr pos  Method[Lewontin|LewontinOnlyHomo|NGSEP]");
				} catch (Exception e) {
					//System.out.println("Try: java -jar ld.jar [VCFldcalculator | 1] [path_vcf] Method[Lewontin|LewontinOnlyHomo|NGSEP]" + e);
				}
			}
			else if (opcion.compareTo("VCFldcalculatorDosage") == 0 || opcion.compareTo("1") == 0) {
				try {
					VCFLDCalculationDosage ld = new VCFLDCalculationDosage();
					
					String ploidy=args[3];
					
					if (ploidy.compareTo("")!=0) {
						ld.VCFldcalulationAll(args[1], args[2], 2);
						ld=null;
					}else if (Integer.parseInt(ploidy)>2){
						ld.VCFldcalulationAll(args[1], args[2], Integer.parseInt(ploidy));
						ld=null;
					}	
					
				} catch (Exception e) {
					System.out.println("Try: java -jar ld.jar [VCFldcalculatorDosage|1] [path_vcf] LewontinOnlyHomoDosage ploidy" + e);
				}
			}
			else if (opcion.compareTo("VCFldcalculator") == 0 || opcion.compareTo("1") == 0) {
				try {
					VCFLDCalculation ld = new VCFLDCalculation();
					ld.VCFldcalulationAll(args[1], args[2]);
					ld=null;
				} catch (Exception e) {
					System.out.println("Try: java -jar ld.jar [VCFldcalculator|1] [path_vcf] Method[Lewontin|LewontinOnlyHomo|NGSEP]" + e);
				}
			}
 
			else if (opcion.compareTo("VCFldcalculatorPairse") == 0 | opcion.compareTo("2") == 0) {
				try {
					VCFLDCalculation ld = new VCFLDCalculation();
					ld.VCFldcalulation(args[1], args[2], args[3], args[4]);
					ld=null;
				} catch (Exception e) {
					System.out.println("Try: java -jar ld.jar [VCFldcalculatorPairse | 2] [path_vcf] Chr pos  Method[Lewontin|LewontinOnlyHomo|NGSEP] " +e);
				}
			}
		
		} catch (Exception e) {
			System.out.println("Try java -jar ld.jar -help");
			System.out.println("Try: java -jar ld.jar [VCFldcalculator | 1] [path_vcf] Method[Lewontin|LewontinOnlyHomo|NGSEP]" + e);
			System.out.println("Try: java -jar ld.jar [ldLewontin | 2] [path_vcf] Chr pos  Method[Lewontin|LewontinOnlyHomo|NGSEP] " +e);
			
		}
		
	}

}
