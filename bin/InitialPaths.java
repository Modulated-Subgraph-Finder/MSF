import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.cli.Options;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;

public class InitialPaths {
	static List<GenesInfo> genesInfoList = new ArrayList<GenesInfo>();
	static List<GenesInteractions> genesInteractionsListComplete = new ArrayList<GenesInteractions>();
	public static List<ArrayList<String>> Allpaths = new ArrayList<ArrayList<String>>();
	static List<String> uniqueGenesInInteraction = new ArrayList<String>();
	static String Step1Output;
	static String network_File;
	static Map<String, Double> genesKeyPValuesPair = new HashMap<String, Double>();
	static Map<String, List<String>> genesKeyInteractionPair = new HashMap<String, List<String>>();
	static GenesInfo geneInfo = new GenesInfo();
	static long threshold = 2L;
	static int rowCounter = 0;
	static int extensionLimit;
	static int mergeLimit;
	static Boolean printExtraFiles = false;
	static String initialPathOutputFile, pvalueInput, interctionInput, dEgType, extraFiles;

	public static void main(String[] args) throws Exception {
		Options options = new Options();
		options.addOption("p", true, "DEG input file path");
		options.addOption("i", true, "Interaction input file path");
		options.addOption("t", true, "DGEA type DEseq2 or EdgeR");
		options.addOption("e", true, "Extension Limit to be used");
		options.addOption("m", true, "Merging limit to be used");
		options.addOption("o", true, "Output folder path");
		options.addOption("k", true, "Extra output files printed");

		CommandLineParser parser = new DefaultParser();
		org.apache.commons.cli.CommandLine cmd = parser.parse(options, args);

		if (cmd.hasOption("m")) {
			mergeLimit = Integer.parseInt(cmd.getOptionValue("m"));
		} else {
			mergeLimit = 1;
		}
		if (cmd.hasOption("p")) {
			pvalueInput = cmd.getOptionValue("p");
		} else {
			System.out.println("Please give DEG input file path ");
			return;
		}
		if (cmd.hasOption("i")) {
			interctionInput = cmd.getOptionValue("i");
		} else {
			System.out.println("Please give Interaction input file path ");
			return;
		}

		if (cmd.hasOption("o")) {
			initialPathOutputFile = cmd.getOptionValue("o");
			network_File = cmd.getOptionValue("o");
		} else {
			System.out.println("Please give Output folder path");
			return;
		}

		if (cmd.hasOption("e")) {
			extensionLimit = Integer.parseInt(cmd.getOptionValue("e"));
		} else {
			extensionLimit = 2;
		}

		if (cmd.hasOption("k")) {
			extraFiles = cmd.getOptionValue("k");
			if (extraFiles.equalsIgnoreCase("yes")) {
				printExtraFiles = true;
			} else {
				printExtraFiles = false;
			}
		}

		if (cmd.hasOption("t")) {
			dEgType = cmd.getOptionValue("t");
		} else {
			System.out.println(
					"Please specify Differential gene expression analysis type DEseq2 or EdgeR with -t option");
			return;
		}
		if ((!dEgType.equalsIgnoreCase("DEseq2")) && (!dEgType.equalsIgnoreCase("EdgeR"))) {
			System.out.println("Please specify Differential gene expression analysis type DEseq2 or EdgeR with -t option");
			return;
		}
		System.out.println("MSF Runing");
		Step1Output = initialPathOutputFile + "InitialGraphs.text";
		network_File = initialPathOutputFile + "NetworkFile.text";
		DataStore.DataStor(pvalueInput, interctionInput, initialPathOutputFile, dEgType);
		geneInfo = DataStore.getpSheet();
		genesInfoList = DataStore.getpSheetList();
		genesInteractionsListComplete = DataStore.getedgeClassListComplete();
		uniqueGenesInInteraction = DataStore.getedgeListSet();
		genesKeyPValuesPair = DataStore.getpsheetKeyValyePair();
		genesKeyInteractionPair = DataStore.getedgeListHashMap();
		identifyingInitialPaths();
		System.out.println("Found Initial Graphs");
		ExtensionMerging.main(Allpaths, initialPathOutputFile, extensionLimit, mergeLimit, printExtraFiles);
	}

	private static List<String> sortpSheetList(List<String> mainValueList) {
		if (mainValueList == null)
			return null;
		List<GenesInfo> psheetList = new ArrayList<GenesInfo>();
		List<String> returnList = new ArrayList<>();
		for (String item : mainValueList) {
			GenesInfo temp = new GenesInfo();
			temp.setGeneName(item);
			try {
				temp.setpValue(genesKeyPValuesPair.get(item));
			} catch (Exception e) {
			}
			if (temp.getpValue() > 0 && temp.getpValue() < 1)
				psheetList.add(temp);
		}
		Collections.sort(psheetList);
		for (GenesInfo returnL : psheetList) {
			returnList.add(returnL.getGeneName());
		}
		return returnList;
	}

	static int returnCounter(int counter, int newCounter) {
		if (newCounter <= 0) {
			return 0;
		} else
			return counter++;
	}

	private static void identifyingInitialPaths() {

		Map<String, List<TempModel>> interactingGeneFrom = new HashMap<String, List<TempModel>>();
		Map<String, List<TempModel>> interactingGeneTo = new HashMap<String, List<TempModel>>();
		HashMap<String, String> geneInteractingAncestor = new HashMap<>();
		HashMap<String, List<String>> geneInteractingAncestorTemp = new HashMap<>();
		DataStore.createHashListFromEdegeClassListTemp(0);
		interactingGeneFrom = DataStore.getInteractingGeneFrom();
		interactingGeneTo = DataStore.getInteractingGeneTo();
		Collections.sort(genesInfoList);
		Set<String> iterationSet = new HashSet<String>(); // all checked genes
		for (GenesInfo psheet : genesInfoList) {
			String rootGeneName = "";
			double combinePValue = 0.0;
			ArrayList<String> selectedGenes = new ArrayList<String>(); // parents
																		// genes
			rootGeneName = psheet.getGeneName();
			geneInteractingAncestor.put(rootGeneName, "");
			if (iterationSet.contains(rootGeneName))
				continue;
			combinePValue = genesKeyPValuesPair.get(rootGeneName);
			iterationSet.add(rootGeneName);
			if (combinePValue >= 0.99)
				continue;
			List<String> interactingGenes = genesKeyInteractionPair.get(rootGeneName);
			if (interactingGenes == null || interactingGenes.contains(null)) {
				int breakpsoint = 1;
				breakpsoint = 2;
				continue;
			}
			for (String string : interactingGenes) {
				geneInteractingAncestor.put(string, rootGeneName);
			}
			interactingGenes = sortpSheetList(interactingGenes);
			selectedGenes.add(rootGeneName);
			int counterCheck = 0;
			for (int counter = counterCheck; interactingGenes != null
					&& counter < interactingGenes.size(); counter = returnCounter(counter, counterCheck)) {
				String genekey = interactingGenes.get(counter);
				interactingGenes.remove(counter);
				counterCheck++;
				if (iterationSet.contains(genekey))
					continue;
				double genePValue = genesKeyPValuesPair.get(genekey);
				if (genePValue >= 1.0)
					continue;
				List<Double> pavlueListString = new ArrayList<>();
				pavlueListString.clear();
				for (String string : selectedGenes) {
					pavlueListString.add(genesKeyPValuesPair.get(string));
				}
				pavlueListString.add(genePValue);
				double calculationResult = GenericFunctions.hartungFunction(pavlueListString);
				if (calculationResult >= combinePValue)
					continue;
				iterationSet.add(genekey);
				combinePValue = calculationResult;
				selectedGenes.add(genekey);
				List<String> newStringListTemp = genesKeyInteractionPair.get(genekey);
				List<String> newStringList = new ArrayList<String>();
				for (String temp : newStringListTemp) {
					if (iterationSet.contains(temp))
						continue;
					newStringList.add(temp);
					if (!geneInteractingAncestor.containsKey(temp))
						geneInteractingAncestor.put(temp.trim(), genekey.trim());
					List<String> tempStringList = geneInteractingAncestorTemp.get(temp.trim());
					if (tempStringList == null || tempStringList.size() == 0)
						tempStringList = new ArrayList<>();
					tempStringList.add(genekey.trim());
					geneInteractingAncestorTemp.put(temp.trim(), tempStringList);
				}
				for (String interact : interactingGenes) {
					if (iterationSet.contains(interact))
						continue;
					newStringList.add(interact);
				}
				interactingGenes = sortpSheetList(newStringList);
				counterCheck = 0;
			}
			if (selectedGenes.size() > 2) {
				for (String stringX : selectedGenes) {
					try {
						String parentGeneOFChildGene = geneInteractingAncestor.get(stringX);
						if (!selectedGenes.contains(parentGeneOFChildGene)) {
							List<String> tempStringList = geneInteractingAncestorTemp.get(stringX);
							for (String string : tempStringList) {
								if (selectedGenes.contains(string)) {
									parentGeneOFChildGene = string;
									break;
								}
							}
						}
						List<TempModel> tempModelList = interactingGeneTo.get(parentGeneOFChildGene);
						if (tempModelList == null)
							continue;
						List<TempModel> xx = tempModelList.stream()
								.filter(x -> x.interactingGene.equalsIgnoreCase(stringX.trim()))
								.collect(Collectors.toList());
						TempModel tempp = new TempModel();
						ArrayList<String> printa = new ArrayList<>();
						String printx = "";
						if (xx.size() < 1) {
							tempModelList = interactingGeneFrom.get(parentGeneOFChildGene);
							xx = tempModelList.stream().filter(x -> x.interactingGene.equalsIgnoreCase(stringX.trim()))
									.collect(Collectors.toList());
							if (xx.size() > 0) {
								tempp = xx.get(0);
								printx = tempp.getInteractingGene() + " " + parentGeneOFChildGene + " "
										+ tempp.getSymbol();
								if (tempp.getSymbol().equals(CustomEnum.activation))
									printx = printx + " 1 " + " 0 ";
								else if (tempp.getSymbol().equals(CustomEnum.Inhibition))
									printx = printx + " 1 " + " 0 ";
								else if (tempp.getSymbol().equals(CustomEnum.reverseInhibitor))
									printx = printx + " 0 " + " 1 ";
								else if (tempp.getSymbol().equals(CustomEnum.reverseActivation))
									printx = printx + " 0 " + " 1 ";
								else if (tempp.getSymbol().equals(CustomEnum.inhibitionActivation)
										|| tempp.getSymbol().equals(CustomEnum.inhibitionInhibition)
										|| tempp.getSymbol().equals(CustomEnum.activationInhibition)
										|| tempp.getSymbol().equals(CustomEnum.activationActivation))
									printx = printx + " 1 " + " 1 ";
							} else {
					
							}

						} else {
							tempp = xx.get(0);
							printx = parentGeneOFChildGene + " " + tempp.getInteractingGene() + " " + tempp.getSymbol();
							if (tempp.getSymbol().equals(CustomEnum.activation))
								printx = printx + " 1 " + " 0 ";
							else if (tempp.getSymbol().equals(CustomEnum.Inhibition))
								printx = printx + " 1 " + " 0 ";
							else if (tempp.getSymbol().equals(CustomEnum.reverseInhibitor))
								printx = printx + " 0 " + " 1 ";
							else if (tempp.getSymbol().equals(CustomEnum.reverseActivation))
								printx = printx + " 0 " + " 1 ";
							else if (tempp.getSymbol().equals(CustomEnum.inhibitionActivation)
									|| tempp.getSymbol().equals(CustomEnum.inhibitionInhibition)
									|| tempp.getSymbol().equals(CustomEnum.activationInhibition)
									|| tempp.getSymbol().equals(CustomEnum.activationActivation))
								printx = printx + " 1 " + " 1 ";
						}
						writeInFIle_NetworkFile(printx);

					} catch (Exception e) {
					}

				}
				Allpaths.add(selectedGenes);
				if (printExtraFiles == true) {
					writeInFIle(selectedGenes, combinePValue);
				}
			}
		}
	}

	private static void writeInFIle_NetworkFile(String output) {

		BufferedWriter bw = null;
		try {
			bw = new BufferedWriter(new FileWriter(network_File, true));
			bw.write(output + "\n");
			bw.newLine();
			bw.flush();
		} catch (IOException ioe) {
			ioe.printStackTrace();
		} finally { // always close the file
			if (bw != null)
				try {
					bw.close();
				} catch (IOException ioe2) {
					// just ignore it
				}
		} // end try/c

	}

	private static synchronized void writeInFIle(ArrayList<String> output, double combineValue) {
		Charset charset = StandardCharsets.UTF_8;
		try {
			Files.write(Paths.get(Step1Output), (output.toString() + " " + combineValue + "\n").getBytes(charset),
					StandardOpenOption.CREATE, StandardOpenOption.APPEND);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}