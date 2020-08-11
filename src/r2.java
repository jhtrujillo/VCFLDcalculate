import java.util.ArrayList;
import java.util.HashMap;

public class r2 {

	// ArrayList<String> haplotypes = new ArrayList<String>();
	HashMap<String, Integer> haplotypes = new HashMap<String, Integer>();
	int numIndHaplotipos = 0;
	double p1 = 0.0;
	double q1 = 0.0;

	public void gethaplotypes(String vcfpath, String SNP1, String SNP2) {

		// --------------------------------------------------------------------------------//
		// Defino las variables a utilizar.
		// --------------------------------------------------------------------------------//
		archivos ar = new archivos();
		String[] datos = ar.leerfichero2(vcfpath);

		String snp1 = "";
		String snp2 = "";

		String A = "";
		String a = "";
		String B = "";
		String b = "";

		String h1 = "";
		String h2 = "";
		String h3 = "";
		String h4 = "";

		int posSNP1 = 0;
		int posSNP2 = 0;

		// --------------------------------------------------------------------------------//
		// Me salto las lineas que tienen comentario.
		// --------------------------------------------------------------------------------//
		boolean escomentario = true;
		int cont = 0;
		while (escomentario) {
			if (datos[cont].toString().contains("#")) {
				cont++;
			} else {
				escomentario = false;
			}
		}

		// --------------------------------------------------------------------------------//
		// Hubico los snps de interes y hago los calculos requiridos para R2 del LD.
		// --------------------------------------------------------------------------------//
		for (int i = cont; i < ar.numerolineas; i++) {
			String[] split = datos[i].split("	");
			int numInd = split.length;
			String snp = split[1];
			String chr = split[0];
			String ref = split[3];
			String alt = split[4];
			double minMaf = Double.parseDouble(split[7].split(";")[2].replace("MAF=", ""));

			// --------------------------------------------------------------------------------//
			// caclulo de los alelos del primer snps
			// --------------------------------------------------------------------------------//
			if (SNP1.compareTo(chr + "_" + snp) == 0) {
				snp1 = chr + " " + snp + " " + ref + " " + alt;
				A = ref;
				a = alt;
				posSNP1 = i;
			}

			// --------------------------------------------------------------------------------//
			// caclulo de los alelos del segundo snps
			// --------------------------------------------------------------------------------//
			if (SNP2.compareTo(chr + "_" + snp) == 0) {
				snp2 = chr + " " + snp + " " + ref + " " + alt;
				B = ref;
				b = alt;
				posSNP2 = i;
				// i = ar.numerolineas;
			}
		}

		// System.out.println(snp1);
		// System.out.println(snp2);

		// --------------------------------------------------------------------------------//
		// Agrego los posibles haplotipos que se pueden generar.
		// --------------------------------------------------------------------------------//
		haplotypes.put(A + "" + B, 0);
		haplotypes.put(a + "" + b, 0);
		haplotypes.put(A + "" + b, 0);
		haplotypes.put(a + "" + B, 0);

		h1 = A + "" + B;
		h2 = a + "" + b;
		h3 = A + "" + b;
		h4 = a + "" + B;

		// --------------------------------------------------------------------------------//
		// ahora realizo el conteo de los haplotipos presentes en la población.
		// --------------------------------------------------------------------------------//
		int numInd = datos[posSNP1].split("	").length - 9;

		// System.out.println(numInd);

		for (int i = 9; i < numInd + 9; i++) {
			String gen1 = datos[posSNP1].split("	")[i].split(":")[0].replace("0", A).replace("1", a);
			String gen2 = datos[posSNP2].split("	")[i].split(":")[0].replace("0", B).replace("1", b);

			String hA = gen1.split("/")[0];
			String ha = gen1.split("/")[1];

			String hB = gen2.split("/")[0];
			String hb = gen2.split("/")[1];

			// --------------------------------------------------------------------------------//
			// Evaluo si ambos SNPs son homocigotos para el individuo.
			// En cualquier otro caso, ya sean valores perdidos o SNPs con ambos
			// heterocigotos no se cuentan.
			// --------------------------------------------------------------------------------//

			if (hA.compareTo(".") != 0 & hB.compareTo(".") != 0) {
				
				if (hA.compareTo(ha) == 0) {
					p1++;
					// Sumo dos por que se generan dos copias en el haplotipo, genoma diploide.
					if (hB.compareTo(hb) == 0) {
						String haplotipo = hA + hB;
						haplotypes.put(haplotipo, haplotypes.get(haplotipo) + 2);
						numIndHaplotipos++;
					}

				}
			}

			

			// --------------------------------------------------------------------------------//
			// Evaluo si tenemos un homocigoto y otro heterocigoto,
			// SNP 1 es homocigoto, pero SNP dos es heterocigoto.
			// SNP 1 es heterocigoto, pero SNP dos es homocigoto.
			// --------------------------------------------------------------------------------//
			if ((hA.compareTo(ha) == 0 && hB.compareTo(hb) != 0 || hA.compareTo(ha) != 0 && hB.compareTo(hb) == 0)
					& hA.compareTo(".") != 0 & hB.compareTo(".") != 0) {

				String haplotipo1 = hA + hB;
				String haplotipo2 = ha + hb;

				haplotypes.put(haplotipo1, haplotypes.get(haplotipo1) + 1);
				haplotypes.put(haplotipo2, haplotypes.get(haplotipo2) + 1);

				numIndHaplotipos++;

				if (hA.compareTo(ha) == 0)
					p1++;
				if (hB.compareTo(hb) == 0)
					q1++;

			}

		}

		System.out.println(p1 + " " + q1 + " " + numIndHaplotipos + "");

		// --------------------------------------------------------------------------------//
		// Calculo la frecuencia de cada uno de los haplotipos
		// --------------------------------------------------------------------------------//
		Double frecuency_h1 = Double.valueOf(haplotypes.get(h1)) / Double.valueOf(numIndHaplotipos * 2);
		Double frecuency_h2 = Double.valueOf(haplotypes.get(h2)) / Double.valueOf(numIndHaplotipos * 2);
		Double frecuency_h3 = Double.valueOf(haplotypes.get(h3)) / Double.valueOf(numIndHaplotipos * 2);
		Double frecuency_h4 = Double.valueOf(haplotypes.get(h4)) / Double.valueOf(numIndHaplotipos * 2);

		// --------------------------------------------------------------------------------//
		// Calculo el valor D, y el valor R2 que representa en LD
		// --------------------------------------------------------------------------------//
		double D = frecuency_h1 - p1 * q1;
		double R2 = D * D / p1 * q1 * (1 - p1) * (1 - q1);

		double p00 = frecuency_h1;
		double p01 = p1;
		double p02 = q1;
		double d = frecuency_h1 - p1 * q1;
		double dPrime;
		if (p01 == 0 || p02 == 0 || p01 == 1 || p02 == 1)
			dPrime = 0;
		else if (d < 0)
			dPrime = d / Math.min(p01 * p02, (1 - p01) * (1 - p02));
		else
			dPrime = d / Math.min(p01 * (1 - p02), (1 - p01) * p02);
		double r2 = d * d;
		if (p01 == 0 || p02 == 0 || p01 == 1 || p02 == 1)
			r2 = 0;
		else
			r2 /= (p01 * p02 * (1 - p01) * (1 - p02));

		R2 = r2;

		// --------------------------------------------------------------------------------//
		// Imprimo los resutlados finales
		// --------------------------------------------------------------------------------//
		
		  System.out.println("Num_Ind\t\t\t" + numInd); System.out.println("SNPs\t\t\t"
		  + snp1 + "," + snp2); System.out.println("Frecuency p1\t\t" + p1);
		  System.out.println("Frecuency q1\t\t" + q1);
		  System.out.println("Num_Ind_counted\t\t" + numIndHaplotipos);
		  System.out.println("Frecuency haplotypes\t" + h1 + ":" + frecuency_h1);
		  System.out.println("\t\t\t" + h2 + ":" + frecuency_h2);
		  System.out.println("\t\t\t" + h3 + ":" + frecuency_h3);
		  System.out.println("\t\t\t" + h4 + ":" + frecuency_h4);
		  System.out.println("D\t\t\t" + D); System.out.println("R2\t\t\t" + R2);
		  System.out.println("");
		 

		System.out.println(SNP1 + "\t" + SNP2 + "\t" + R2 + "\t" + D);

	}

	public void ld_multiple(String vcfpath, String SNP1) {
		archivos ar = new archivos();
		String[] datos = ar.leerfichero2(vcfpath);

		boolean escomentario = true;
		int cont = 0;
		while (escomentario) {
			if (datos[cont].toString().contains("#")) {
				cont++;
			} else {
				escomentario = false;
			}
		}

		for (int i = cont; i < ar.numerolineas; i++) {
			String[] split = datos[i].split("	");
			String SNP2 = split[0] + "_" + split[1];
			r2 r2ld = new r2();
			r2ld.gethaplotypes(vcfpath, SNP1, SNP2);
			r2ld = null;
		}

		ar = null;
	}

	public static void main(String[] args) {
		r2 r2ld = new r2();

		try {
			String opcion = args[0];

			if (opcion.compareTo("-h") == 0) {
				System.out.println("Try java -jar ld.jar 1 vcf SNP1 SNP2");
				System.out.println("Try java -jar ld.jar 2 vcf SNP");
			}

			else if (opcion.compareTo("1") == 0) {
				try {
					r2ld.gethaplotypes(args[1], args[2], args[3]);
				} catch (Exception e) {
					System.out.println("Try java -jar ld.jar 1 vcf SNP1 SNP2");
				}
			}

			else if (opcion.compareTo("2") == 0) {
				try {
					r2ld.ld_multiple(args[1], args[2]);
				} catch (Exception e) {
					System.out.println("Try java -jar ld.jar 2 vcf SNP");
				}
			}

		} catch (Exception e) {
			System.out.println("Try java -jar fingerprint.jar -help");
			System.out.println("Try java -jar ld.jar 1 vcf SNP1 SNP2");
			System.out.println("Try java -jar ld.jar 2 vcf SNP");
		}

		// r2ld.ld_multiple("/home/estuvar4/Desktop/mergevcf.95ids.b_fourthfiltered.vcf",
		// "scaffold_16195_79974");

		r2ld.gethaplotypes("/home/estuvar4/Desktop/mergevcf.95ids.b_fourthfiltered.vcf", "scaffold_12997_216675",
				"scaffold_12997_217253");

		/*
		 * 
		 * r2ld.gethaplotypes(
		 * "/home/estuvar4/Desktop/mergevcf.95ids.b_fourthfiltered.vcf",
		 * "scaffold_7561_7519", "scaffold_12997_216675");
		 * 
		 * r2ld.gethaplotypes(
		 * "/home/estuvar4/Desktop/mergevcf.95ids.b_fourthfiltered.vcf",
		 * "scaffold_16195_79974", "scaffold_16195_80046");
		 * 
		 * r2ld.gethaplotypes(
		 * "/home/estuvar4/Desktop/mergevcf.95ids.b_fourthfiltered.vcf",
		 * "scaffold_16195_60995", "scaffold_16195_63677");
		 * 
		 * r2ld.gethaplotypes(
		 * "/home/estuvar4/Desktop/mergevcf.95ids.b_fourthfiltered.vcf",
		 * "scaffold_12298_420354", "scaffold_12298_451660");
		 * 
		 * r2ld.gethaplotypes(
		 * "/home/estuvar4/Desktop/mergevcf.95ids.b_fourthfiltered.vcf",
		 * "scaffold_5699_69673", "scaffold_5699_114979");
		 */

	}
}
