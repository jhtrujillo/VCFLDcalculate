import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import ngsep.variants.CalledGenomicVariant;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFRecord;

public class VCFldcalculator {

	// --------------------------------------------
	// Cargo el VCF a memoria con modulos de NGSEP
	// --------------------------------------------
	public List<VCFRecord> loadVCF(String filename, String SNP1, String SNP2) throws IOException {
		VCFFileReader in = new VCFFileReader(filename);
		List<VCFRecord> recordsInMemory = new LinkedList<>();

		VCFRecord record1 = null;
		VCFRecord record2 = null;

		// ----------------------------------------------------
		// Busco los SNPs y los alamanceno en dos VCFrecords.
		// ----------------------------------------------------
		Iterator<VCFRecord> it = in.iterator();
		while (it.hasNext()) {
			VCFRecord record = it.next();
			List<CalledGenomicVariant> calls = record.getCalls();

			String snp = calls.get(0).getSequenceName() + "_" + calls.get(0).getFirst();

			if (snp.compareTo(SNP1) == 0) {
				record1 = record;
				recordsInMemory.add(record1);
				// System.out.println(calls.get(0).getSequenceName() + " " +
				// calls.get(0).getFirst());
			}

			if (snp.compareTo(SNP2) == 0) {
				record2 = record;
				recordsInMemory.add(record2);
				break;
				// System.out.println(calls.get(0).getSequenceName() + " " +
				// calls.get(0).getFirst());
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

		double p1 = 0.0;
		double q1 = 0.0;
		double freHap = 0.0;
		int numInd = 0;
		int cont_hapSNP = 0;

		// --------------------------------------------
		// Frecuencia de los alelos menores y de haplotipo AB.
		// --------------------------------------------
		for (int i = 0; i < n; i++) {
			CalledGenomicVariant call1 = calls1.get(i);
			CalledGenomicVariant call2 = calls2.get(i);

			if (call1.isUndecided() || call2.isUndecided())
				continue;
			if (call1.isHeterozygous() & call2.isHeterozygous())
				continue;

			numInd++;

			if (call1.isHomozygousReference()) {
				p1++;
				if (call2.isHomozygousReference()) {
					cont_hapSNP++;
				}
				else if (call2.isHeterozygous()) {
					String aleloRefSnp = call2.getReference();
					String aleloRefGen = call2.getAlleles()[0];

					if (aleloRefSnp.compareTo(aleloRefGen) == 0) {
						cont_hapSNP++;
					}
				}
			}

			if (call2.isHomozygousReference()) {
				q1++;
			}

		}

		// --------------------------------------------
		// Haplotipo a buscar
		// --------------------------------------------
		String refSNP1 = calls1.get(0).getReference();
		String refSNP2 = calls2.get(0).getReference();
		String hapSNP = refSNP1 + "" + refSNP2;

		// --------------------------------------------
		// Calculos estadisticos para el R2
		// --------------------------------------------
		p1 = p1 / numInd;
		q1 = q1 / numInd;
		double frecHap = cont_hapSNP / numInd;
		double D = frecHap - p1 * q1;
		double R2 = D * D / p1 * q1 * (1 - p1) * (1 - q1);

		Double dPrime;

		if (p1 == 0 || q1 == 0 || p1 == 1 || q1 == 1)
			dPrime = 0.0;
		else if (D < 0)
			dPrime = D / Math.min(p1 * q1, (1 - p1) * (1 - q1));
		else
			dPrime = D / Math.min(p1 * (1 - q1), (1 - p1) * q1);
		double r2 = D * D;
		if (p1 == 0 || q1 == 0 || p1 == 1 || q1 == 1)
			r2 = 0;
		else
			r2 /= (p1 * q1 * (1 - p1) * (1 - q1));

		System.out.print("NumInd\t");
		// System.out.print("p1\t");
		// System.out.print("q1\t");
		// System.out.print("Frec\t");
		System.out.print("D\t");
		System.out.print("Dprime\t");
		System.out.print("R2\t");
		System.out.println("r2\t");

		System.out.print(numInd + "\t");
		// System.out.print(p1+"\t");
		// System.out.print(q1+"\t");
		// System.out.print(p1 * q1+"\t");
		System.out.print(D + "\t");
		System.out.print(dPrime + "\t");
		System.out.print(R2 + "\t");
		System.out.println(r2 + "\t");

	}

	public static void main(String[] args) throws IOException {
		//VCFldcalculator ld = new VCFldcalculator();
		//ld.CalculateLD("/home/estuvar4/Desktop/mergevcf.95ids.b_fourthfiltered.vcf", "scaffold_12997_216675",
		//	"scaffold_12997_217253");

		VCFldcalculator r2ld = new VCFldcalculator();

		try {
			String opcion = args[0];

			if (opcion.compareTo("-h") == 0) {
				System.out.println("Try java -jar ld.jar 1 vcf SNP1 SNP2");
				System.out.println("Try java -jar ld.jar 2 vcf SNP");
			}

			else if (opcion.compareTo("1") == 0) {
				try {
					r2ld.CalculateLD(args[1], args[2], args[3]);
				} catch (Exception e) {
					System.out.println("Try java -jar r2ld.jar 1 vcf SNP1 SNP2");
				}
			}

		} catch (Exception e) {
			System.out.println("Try java -jar fingerprint.jar -help");
			System.out.println("Try java -jar ld.jar 1 vcf SNP1 SNP2");
		}

	}
}
