import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;
import ngsep.variants.VariantCallReport;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFLDCalculator;
import ngsep.vcf.VCFRecord;

public class VCFLDCalculationDosage {
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
	double homcigousLeve=0.3;

	public float roundToArray(float value, float [] array){
		 
		 float best_abs = 1.0f;
		 float rounded = value;
		 
		 for(int i=0; i < array.length; i++){
			 if(best_abs > Math.abs(value-array[i]))
			 {
				 best_abs = value-array[i];
				 rounded = array[i];
			 }		 
		 }
		 return rounded;
	 }
	
	// --------------------------------------------
		// Calcular valores de LD segun On Measures of Gametic Disequilibrium R. C.
		// Lewontin. Solo contando los sitios dobles homocigotos.
		// --------------------------------------------
		public void getValuesLDLewontinOnlyHomo_Dosage(VCFRecord record1, VCFRecord record2, int ploidy) throws IOException {
			
			List<CalledGenomicVariant> calls1 = record1.getCalls();
			List<CalledGenomicVariant> calls2 = record2.getCalls();

			String SNP1 = calls1.get(0).getSequenceName() + " " + calls1.get(0).getFirst();
			String SNP2 = calls2.get(0).getSequenceName() + " " + calls2.get(0).getFirst();
			
			GenomicVariant var1 = record1.getVariant();
			String [] alleles1 = var1.getAlleles();
			
			

			int n = calls1.size();

			p1 = 0.0;
			p2 = 0.0;
			q1 = 0.0;
			q2 = 0.0;

			freHap = 0.0;
			frecAB = 0.0;
			Dprime = 0;
			numInd = 0;
			
			
			if (ploidy<2) {
				ploidy=2;
			}
			
			//System.out.println("Ploidia "+ploidy  );
			
			float ploidyLevels[] = new float[ploidy+1];

			for(int y=0; y <= ploidy;y++){
				ploidyLevels[y] = (1.0f/ploidy) * y;
				//System.out.println(ploidyLevels[y] );
			}
			
			
			// System.out.println("Num Ind Lewtion inicio "+ numInd);

			// System.out.println("CC n "+n);

			// --------------------------------------------
			// Frecuencia de los alelos menores y de haplotipo AB.
			// --------------------------------------------
			for (int i = 0; i < n; i++) {
				
				
				GenomicVariant var2 = record2.getVariant();
				String [] alleles2 = var2.getAlleles();
				
				CalledGenomicVariant call1 = calls1.get(i);
				CalledGenomicVariant call2 = calls2.get(i);
					
					
				if (call1.isUndecided() || call1.isHeterozygous())	continue;
				if (call2.isUndecided() || call2.isHeterozygous()) continue;

				VariantCallReport report1 = call1.getCallReport();
				VariantCallReport report2 = call2.getCallReport();
			
	 			float countRef1 = report1.getCount(alleles1[0]);
	    		float countAlt1 = report1.getCount(alleles1[1]);
	    		
	    		float countRef2 = report2.getCount(alleles2[0]);
	    		float countAlt2 = report2.getCount(alleles2[1]);
				
	    		float dosage1=0;
	    		float dosage2=0;
	    		
	    		//Depends of ploidy assign a value to dosage
	    		if((countRef1 + countAlt1) > 0){
	    			float dosage1_tmp = countRef1 / (countRef1 + countAlt1);
	    			dosage1 = roundToArray(dosage1_tmp, ploidyLevels);
	    			//System.out.println(dosage1_tmp);
    			}
	    		
	    		//Depends of ploidy assign a value to dosage
	    		if((countRef2 + countAlt2) > 0){
	    			float dosage2_tmp = countRef2 / (countRef2 + countAlt2);
	    			dosage2 = roundToArray(dosage2_tmp, ploidyLevels);
    			}
	    		//System.out.println(countRef1+" "+countAlt1+" "+dosage1+" || "+countRef2+" "+countAlt2+" "+dosage2);
	    		//System.out.println(dosage1+" "+dosage2);
	    		
				numInd++;

				// Dobles homocigotos.
				if (dosage1 >= this.homcigousLeve) {
					p1 += 1;
				}

				// Dobles homocigotos.
				if (dosage2 >= this.homcigousLeve) {
					q1 += 1;
				}

				if ( dosage1 >= this.homcigousLeve && dosage2 >= this.homcigousLeve) {
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
			
			if (SNP1.compareTo(SNP2)!=0) {
				System.out.println(SNP1 + "\t" + SNP2 + "\t" + D + "\t" + dPrime + "\t" + r2 + "\t" + numInd+"\t"+ploidy);
			}

			

			this.r2_lewinton = r2;

		}
	
	
	public void CalculateLDNGSEP(VCFRecord record1, VCFRecord record2) throws IOException {
		VCFLDCalculator ld = new VCFLDCalculator();
		ld.calculateLDStatistics(record1, record2);

	}

	public void VCFldcalulation(String filename, String Chr, String pos, String opcion, int ploidy) throws IOException {

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
		System.out.println("SNP1" + "\t" + "SNP2" + "\t" + "D" + "\t" + "dPrime" + "\t" + "r2" + "\t" + "numInd"+ "\t" + "Poidy");

		while (it2.hasNext()) {
			VCFRecord SNPrecord2 = it2.next();

			if (opcion.compareTo("LewontinOnlyHomoDosage") == 0) {
				getValuesLDLewontinOnlyHomo_Dosage(SNPrecord1, SNPrecord2, ploidy);
			}
			

		}

	}

	public void VCFldcalulationAll(String filename, String opcion, int ploidy) throws IOException {

		VCFFileReader in = new VCFFileReader(filename);
		Iterator<VCFRecord> it = in.iterator();

		System.out.println("SNP1" + "\t" + "SNP2" + "\t" + "D" + "\t" + "dPrime" + "\t" + "r2" + "\t" + "numInd"+ "\t" + "Poidy");

		// Saco el VCFRecord del SNP de interes.
		while (it.hasNext()) {
			VCFRecord SNPrecord1 = it.next();

			// Comparo el VCFRecord seleccionado con el resto de SNPs.
			VCFFileReader in2 = new VCFFileReader(filename);
			Iterator<VCFRecord> it2 = in2.iterator();

			while (it2.hasNext()) {
				VCFRecord SNPrecord2 = it2.next();

				 if (opcion.compareTo("LewontinOnlyHomoDosage") == 0) {
					getValuesLDLewontinOnlyHomo_Dosage(SNPrecord1, SNPrecord2,  ploidy);
				} else if (opcion.compareTo("NGSEP") == 0) {
					CalculateLDNGSEP(SNPrecord1, SNPrecord2);
				}
			}
		}
	}

	public static void main(String[] args) throws IOException {

		
		 VCFLDCalculationDosage ld = new VCFLDCalculationDosage();
		 
	
		String Chr = "scaffold_12997";
		String snp = "216675";
		
		
		
		ld.VCFldcalulation("/home/estuvar4/Desktop/tmp.vcf", Chr, snp, "LewontinOnlyHomoDosage", 10);
		//ld.VCFldcalulation("/home/estuvar4/Desktop/tmp.vcf", Chr, snp, "LewontinOnlyHomoDosage", 2);
		
		//ld.VCFldcalulation("/home/estuvar4/Desktop/tmp.vcf", Chr, snp, "LewontinOnlyHomoDosage", 10);
		
		
		//ld.VCFldcalulationAll("/home/estuvar4/Desktop/mergevcf.95ids.b_fourthfiltered.vcf", "LewontinOnlyHomoDosage", ploidy);
	
		
		// ld.VCFldcalulation("/home/estuvar4/Desktop/tmp.vcf", Chr, snp, "Lewontin");

		// ld.VCFldcalulation("/home/estuvar4/Desktop/tmp.vcf", Chr, snp,
		// "LewontinOnlyHomo");

		// ld.VCFldcalulation("/home/estuvar4/Desktop/tmp.vcf", Chr, snp, "NGSEP");
		
		//ld.VCFldcalulationAll("/home/estuvar4/Desktop/tmp.vcf", "Lewontin");
		//ld.VCFldcalulationAll("/home/estuvar4/Desktop/tmp.vcf", "LewontinOnlyHomo");
		//ld.VCFldcalulationAll("/home/estuvar4/Desktop/tmp.vcf", "NGSEP");
		
		
		/*
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
		*/
	}

}
