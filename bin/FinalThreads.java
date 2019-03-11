import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.math3.stat.inference.TTest;


public class FinalThreads {

	static Map<String, List<String>> genesKeyInteractionPair = new HashMap<String, List<String>>();
	static Map<String, List<String>> genesKeyInteractionPairSourceSink = new HashMap<String, List<String>>();
	static Map<String, List<String>> interactingGeneFrom = new HashMap<String, List<String>>();
	static Map<String, List<String>> interactingGeneFromNetwork = new HashMap<String, List<String>>();
	static Map<String, List<String>> interactingGeneTo = new HashMap<String, List<String>>();
	static Map<String, List<String>> interactingGeneToNetwork = new HashMap<String, List<String>>();
	static List<GenesInteractions> genesInteractionsListComplete = new ArrayList<GenesInteractions>();
	static List<String> uniqueGenesInInteraction = new ArrayList<String>();
	static List<GenesInfo> genesInfoList = new ArrayList<GenesInfo>();
	static Map<String, Double> genesKeyPValuesPair = new HashMap<String, Double>();
	static Map<String, String> genesKeyLogFoldChangesPair = new HashMap<String, String>();
	static int rowCounter = 0;
	static Hashtable<String, String> geneInfoHashTable = new Hashtable<>();
	static List<String> genesList = new ArrayList<>();
	static List<String> ReactomeList = new ArrayList<>();
	static String fileoutPathText, fileoutPathTextTemp,sourceWeightagePrintFile;
	static int sourceCounter = 0, sinkCounter = 0;
	static String network_File = "";

	public static void main(List<ArrayList<String>> AllMergedpaths, String network_FileTemp) throws Exception {
		fileoutPathText = DataStore.getOutputPath() + "SourcesAndSinks.text";
		sourceWeightagePrintFile = DataStore.getOutputPath() + "sourceWeightage.text";
		fileoutPathTextTemp = DataStore.getOutputPath() + "Ttest.text";
		network_File = network_FileTemp;
		System.out.println("Finding Sources and Sinks");
		genesInfoList = DataStore.getpSheetList();
		genesInteractionsListComplete = DataStore.getedgeClassListComplete();
		uniqueGenesInInteraction = DataStore.getedgeListSet();
		genesKeyPValuesPair = DataStore.getpsheetKeyValyePair();
		genesKeyLogFoldChangesPair =  DataStore.getpsheetKeyFolderValue();
		genesKeyInteractionPair = DataStore.getedgeListHashMap();
		geneInfoHashTable = DataStore.getgeneInfoHashTable();
		createHashListFromEdegeClassList(0);
		int listNumer = 1;
		for (ArrayList<String> arrayList : AllMergedpaths) {
			if (listNumer != 1) {
				List<String> prepareListTemp = new ArrayList<String>();
			}
			step4(arrayList, (listNumer));
			printSourceWeightage(arrayList);
			listNumer++;
		}
		List<String> prepareListTemp = new ArrayList<String>();
		prepareListTemp.add("SOurce" + sourceCounter);
		prepareListTemp.add("Sink" + sinkCounter);
		System.out.println("Sources and Sinks Identified");
		System.out.println("Creating Output Files");
		System.out.println("DONE");
	}

	private static void printSourceWeightage(ArrayList<String> arrayList) {
		StringBuilder prepareList = new StringBuilder();
		for (String string : arrayList) {
			prepareList.append(string);
			String geneType = checkGeneSourceSink(arrayList, string);
				if (geneType.equalsIgnoreCase("source")) {
					List<String> temp = sourceGeneWeight(arrayList, string);
					Double weightage = (double) (temp.size() / (double) arrayList.size() * 100);
					prepareList.append("  "+String.format("%.3f", weightage*10)+ "  "+genesKeyLogFoldChangesPair.get(string.toLowerCase()));
				}
				else
					prepareList.append("  "+String.format("%.3f", 10.0)+ "  "+genesKeyLogFoldChangesPair.get(string.toLowerCase()));
				prepareList.append("\n");
		}
		Charset charset = StandardCharsets.UTF_8;
		try {
			Files.write(Paths.get(sourceWeightagePrintFile), (prepareList.toString() + "\n").getBytes(charset),
					StandardOpenOption.CREATE, StandardOpenOption.APPEND);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static synchronized void step4(List<String> genesList, int pathNumber) {
		List<String> prepareList = new ArrayList<String>();
		sourceCounter = 0;
		sinkCounter = 0;
		prepareList.add("Graph " + pathNumber);
		writeInFIle(prepareList);
		prepareList = new ArrayList<String>();
		List<Double> pavlueListString = new ArrayList<>();
		for (String string : genesList) {
			ReactomeList.add(string);
			prepareList.add(string);
			pavlueListString.add(genesKeyPValuesPair.get(string));
		}
		writeInFIle(prepareList);
		readNetworkFile(network_File);
		for (String string : genesList) {
			prepareList = new ArrayList<String>();
			prepareList.add(string);
			pavlueListString.clear();
			prepareList.add(geneInfoHashTable.get(string));
			prepareList.add(String.valueOf(genesKeyPValuesPair.get(string)));
			prepareList.add(checkGeneSourceSink(genesList, string));
			String geneType = checkGeneSourceSink(genesList, string);
			{
				if (geneType.equalsIgnoreCase("source")) {
					List<String> temp = sourceGeneWeight(genesList, string);
					Double weightage = (double) (temp.size() / (double) genesList.size() * 100);
					prepareList.add(temp.size() + "/" + genesList.size() + " = " + String.format("%.3f", weightage));
				}
				writeInFIle(prepareList);
			}
			if (geneType.equalsIgnoreCase("source")) {
				List<String> temp = sourceDownSteamGeneWeight(genesList, string);
				List<String> upStreamInteractions = getUpstreamGenesNotinNetworkFile(string);
				List<String> tempUpstreamIneraction = sourceDownSteamGene(genesList, string, upStreamInteractions);
				List<String> tempList = new ArrayList<>();
				System.out.println("gene itself "+string);
				System.out.println("upStreamInteractions size"+ upStreamInteractions.size());
				System.out.println("upStreamInteractions gene"+ upStreamInteractions.toString());
				double[] upstreamArray1 = new double [upStreamInteractions.size()]; 
				int countX= 0;
				for (String string2 : upStreamInteractions) {
					tempList.add(string2);
					upstreamArray1[countX++] = (double) genesKeyPValuesPair.get(string2);
				}
				for (String string2 : tempUpstreamIneraction) { 
					tempList.add(string2);
				}
				if (upStreamInteractions != null && upStreamInteractions.size() > 0) 
					writeInFIleTemp(tempList, tempList.size(), string);
			
				// t test on genes
				List<String> downStreamInteractions = getDownstreamGenesOnlyinNetworkFile(string);				
				List<Double> downsteamPvalues = new ArrayList<>();
				countX= 0;
				System.out.println("downStreamInteractions size"+ downStreamInteractions.size());
				System.out.println("downStreamInteractions gene"+ downStreamInteractions.toString());
				
				double[] downstreamArray1 = new double [downStreamInteractions.size()];
				for (String string2 : downStreamInteractions) {
					downsteamPvalues.add(genesKeyPValuesPair.get(string2));
					downstreamArray1[countX++] = (double) genesKeyPValuesPair.get(string2);
				}
				double p = new TTest().tTest(upstreamArray1,downstreamArray1);
				System.out.println(String.valueOf(p));
				
			} // vlose of source 
		} 
	}

	private static synchronized void writeInFIleTemp(List<String> temp, int size, String gene) {
		Charset charset = StandardCharsets.UTF_8;
		String finalValue = "Gene " + gene + "\n";
		for (int count = 0; count < temp.size() && count < size; count++) {
			finalValue = finalValue + temp.get(count) + " "
					+ String.valueOf(genesKeyPValuesPair.get(temp.get(count)) + " ");
		}
		try {
			Files.write(Paths.get(fileoutPathTextTemp), (finalValue.toString() + "\n").getBytes(charset),
					StandardOpenOption.CREATE, StandardOpenOption.APPEND);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static List<String> sourceGeneWeight(List<String> genesList, String gene) {
		List<String> finalSelectedList = new ArrayList<>();
		finalSelectedList.add(gene);
		HashMap<String, String> whoBroughtWhom = new HashMap<>();
		List<String> interactionGene = getAllInteractionFromList(genesList, gene);
		List<String> doneGnes = new ArrayList<>();
		for (String string : interactionGene)
			whoBroughtWhom.put(string, gene);
		int counterCheck = 0;
		doneGnes.add(gene);
		for (int counter = counterCheck; interactionGene != null
				&& counter < interactionGene.size(); counter = returnCounter(counter, counterCheck)) {
			String genekey = interactionGene.get(counter);
			interactionGene.remove(counter);
			doneGnes.add(genekey);
			counterCheck++;
			String innerGene = checkGeneSourceSink(genesList, genekey);
			if (innerGene.equalsIgnoreCase("sink") || innerGene.equalsIgnoreCase("Intermediate")) {
				if (finalSelectedList.contains(whoBroughtWhom.get(genekey)))
					finalSelectedList.add(genekey);
				if (innerGene.equalsIgnoreCase("Intermediate")) {
					List<String> interactionGeneTemp = getAllInteractionFromList(genesList, genekey);
					for (String string : interactionGeneTemp) {
						if (!doneGnes.contains(string)) {
							interactionGene.add(string);
							whoBroughtWhom.put(string, genekey);
						}
					}
					counterCheck = 0;
				}
			}
		}
		return finalSelectedList;
	}

	private static List<String> sourceDownSteamGeneWeight(List<String> genesList, String gene) {
		List<String> finalSelectedList = new ArrayList<>();
		List<String> finalIteratedList = new ArrayList<>();
		finalSelectedList.add(gene);
		finalIteratedList.add(gene);
		HashMap<String, String> whoBroughtWhom = new HashMap<>();
		List<String> interactionGene = getAllInteractionFromList(genesList, gene);
		List<String> doneGnes = new ArrayList<>();
		for (String string : interactionGene)
			whoBroughtWhom.put(string, gene);
		int counterCheck = 0;
		doneGnes.add(gene);
		for (int counter = counterCheck; interactionGene != null
				&& counter < interactionGene.size(); counter = returnCounter(counter, counterCheck)) {
			String genekey = interactionGene.get(counter);
			interactionGene.remove(counter);
			doneGnes.add(genekey);
			counterCheck++;
			String innerGene = checkGeneSourceSink(genesList, genekey);
			if (innerGene.equalsIgnoreCase("sink") || innerGene.equalsIgnoreCase("Intermediate")) {
				if (finalSelectedList.contains(whoBroughtWhom.get(genekey)))
					finalSelectedList.add(genekey);
				if (innerGene.equalsIgnoreCase("Intermediate")) {
					List<String> interactionGeneTemp = getAllInteractionFromList(genesList, genekey);
					for (String string : interactionGeneTemp) {
						if (!doneGnes.contains(string)) {
							interactionGene.add(string);
							whoBroughtWhom.put(string, genekey);
						}
					}
					counterCheck = 0;
				}
			}
		}
		return finalSelectedList;
	}

	/* old function previously used but no more in use */
	private static List<ArrayList<String>> sourceDownSteamGeneWeightTemp(List<String> genesList, String gene) {
		ArrayList<String> finalSelectedList = new ArrayList<>();
		ArrayList<String> finalSelectedListTemp = new ArrayList<>();
		List<ArrayList<String>> allConnectedPath = new ArrayList<ArrayList<String>>();
		finalSelectedListTemp.add(gene);
		HashMap<String, String> whoBroughtWhom = new HashMap<>();
		List<String> interactionGene = getAllInteractionFromList(genesList, gene);
		List<String> doneGnes = new ArrayList<>();
		for (String string : interactionGene)
			whoBroughtWhom.put(string, gene);
		int counterCheck = 0;
		doneGnes.add(gene);
		for (int counter = counterCheck; interactionGene != null && counter < interactionGene.size(); counter++) {
			String genekey = interactionGene.get(counter);
			interactionGene.remove(counter);
			doneGnes.add(genekey);
			counterCheck++;
			String innerGene = checkGeneSourceSink(genesList, genekey);
			if (innerGene.equalsIgnoreCase("sink") || innerGene.equalsIgnoreCase("Intermediate")) {
				if (finalSelectedListTemp.contains(whoBroughtWhom.get(genekey))) {
					finalSelectedListTemp.add(genekey);
				}
				if (innerGene.equalsIgnoreCase("Intermediate")) {
					List<String> interactionGeneTemp = getAllInteractionFromList(genesList, genekey);
					for (String string : interactionGeneTemp) {
						if (!doneGnes.contains(string)) {
							String innerGene2 = checkGeneSourceSink(genesList, string);
							if (innerGene2.equalsIgnoreCase("sink") || innerGene2.equalsIgnoreCase("Intermediate")) {
								// if (finalSelectedListTemp.contains(whoBroughtWhom.get(string)))
								{
									finalSelectedList.add(finalSelectedListTemp.get(0));
									finalSelectedList.add(finalSelectedListTemp.get(1));
									finalSelectedList.add(string);
									allConnectedPath.add(finalSelectedList);
									finalSelectedList = new ArrayList<>();
								}
							}
						}
					}
				}
			}
		}
		return allConnectedPath;
	}

	private static List<String> sourceDownSteamGene(List<String> genesList, String gene,
			List<String> alreadyProcessed) {
		ArrayList<String> finalSelectedList = new ArrayList<>();
		finalSelectedList.add(gene);
		HashMap<String, String> whoBroughtWhom = new HashMap<>();
		List<String> interactionGene = getDownstreamGenesOnlyinNetworkFile(gene);
		List<String> doneGnes = new ArrayList<>();
		for (String string : interactionGene)
			whoBroughtWhom.put(string, gene);
		int counterCheck = 0;
		doneGnes.add(gene);
		for (int counter = counterCheck; interactionGene != null && counter < interactionGene.size(); counter++) {
			String genekey = interactionGene.get(counter);
			interactionGene.remove(counter);
			counterCheck++;
			if (!alreadyProcessed.contains(genekey) && !doneGnes.contains(genekey)) {
				doneGnes.add(genekey);
				finalSelectedList.add(genekey);
			}
			List<String> interactionGeneTemp = getDownstreamGenesOnlyinNetworkFile(genekey);
			for (String string : interactionGeneTemp) {
				if (!alreadyProcessed.contains(string) && !doneGnes.contains(genekey)) {
					finalSelectedList.add(string);
					doneGnes.add(gene);
				}
			}
		}
		return finalSelectedList;
	}

	private static List<ArrayList<String>> sourceUpstreamSteamGeneWeightTemp(List<String> interactionGenes,
			String gene) {
		ArrayList<String> finalSelectedListTemp = new ArrayList<>();
		List<ArrayList<String>> allConnectedPath = new ArrayList<ArrayList<String>>();
		// finalSelectedListTemp.add(gene);
		for (String string : interactionGenes) {
			finalSelectedListTemp.add(string);
			List<String> nextInteractions = getUpstreamGenesNotinNetworkFile(string);
			for (String string2 : nextInteractions) {
				finalSelectedListTemp.add(string2);
				List<String> tempList = new ArrayList<>(finalSelectedListTemp);
				Collections.reverse(tempList);
				allConnectedPath.add(new ArrayList<>(tempList));
				finalSelectedListTemp = new ArrayList<>();
			}
		}
		return allConnectedPath;
	}

	private static List<ArrayList<String>> sourceUpsteamGeneWeight(String geneItself, List<String> downStreamPAth) {
		ArrayList<String> finalSelectedList = new ArrayList<>();
		ArrayList<String> finalSelectedListTemp = new ArrayList<>();
		List<ArrayList<String>> allConnectedPath = new ArrayList<ArrayList<String>>();
		finalSelectedListTemp.add(geneItself);
		HashMap<String, String> whoBroughtWhom = new HashMap<>();
		List<String> interactionGene = getAllInteractionFromList(genesList, geneItself);
		List<String> doneGnes = new ArrayList<>();
		for (String string : interactionGene)
			whoBroughtWhom.put(string, geneItself);
		int counterCheck = 0;
		doneGnes.add(geneItself);
		for (int counter = counterCheck; interactionGene != null && counter < interactionGene.size(); counter++) {
			String genekey = interactionGene.get(counter);
			interactionGene.remove(counter);
			doneGnes.add(genekey);
			counterCheck++;
			String innerGene = checkGeneSourceSink(genesList, genekey);
			if (innerGene.equalsIgnoreCase("sink") || innerGene.equalsIgnoreCase("Intermediate")) {
				if (finalSelectedListTemp.contains(whoBroughtWhom.get(genekey))) {
					finalSelectedListTemp.add(genekey);
				}
				if (innerGene.equalsIgnoreCase("Intermediate")) {
					List<String> interactionGeneTemp = getAllInteractionFromList(genesList, genekey);
					for (String string : interactionGeneTemp) {
						if (!doneGnes.contains(string)) {
							String innerGene2 = checkGeneSourceSink(genesList, string);
							if (innerGene2.equalsIgnoreCase("sink") || innerGene2.equalsIgnoreCase("Intermediate")) {
								finalSelectedList.add(finalSelectedListTemp.get(0));
								finalSelectedList.add(finalSelectedListTemp.get(1));
								finalSelectedList.add(string);
								allConnectedPath.add(finalSelectedList);
								finalSelectedList = new ArrayList<>();
							}
						}
					}
				}
			}
		}
		return allConnectedPath;
	}

	static int returnCounter(int counter, int newCounter) {
		if (newCounter <= 0) {
			return 0;
		} else
			return counter++;
	}

	private static List<String> getAllInteractionFromList(List<String> genesList, String gene) {
		List<String> allInteractionList = new ArrayList<>();
		GenesInteractions edgeClass = new GenesInteractions();
		edgeClass.getColmn2();
		List<String> sinkInteraction = new ArrayList<>();
		List<String> interctionTempTotally = genesKeyInteractionPairSourceSink.get(gene);
		if (null == interctionTempTotally)
			interctionTempTotally = new ArrayList<>();
		sinkInteraction = interactingGeneFrom.get(gene);
		if (null != sinkInteraction)
			if (sinkInteraction.contains(gene))
				sinkInteraction.remove(gene);
		List<String> sourceIneraction = new ArrayList<>();
		sourceIneraction = interactingGeneTo.get(gene);
		if (null != sourceIneraction)
			if (sourceIneraction.contains(gene))
				sourceIneraction.remove(gene);
		Boolean sink = false, source = false;
		String sinkGeneName = null, sourceGeneName = null;
		if (null != sinkInteraction)
			for (String sinkGene : sinkInteraction) {
				if (genesList.contains(sinkGene) && interctionTempTotally.contains(sinkGene)) {
					allInteractionList.add(sinkGene);
				}
			}
		if (null != sourceIneraction)
			for (String sourceGene : sourceIneraction) {
				if (genesList.contains(sourceGene) && interctionTempTotally.contains(sourceGene)) {
					allInteractionList.add(sourceGene);
				}
			}
		return allInteractionList;
	}

	private static String checkGeneSourceSink(List<String> genesList, String gene) {
		GenesInteractions edgeClass = new GenesInteractions();
		edgeClass.getColmn2();
		List<String> sinkInteraction = new ArrayList<>();
		List<String> interctionTempTotally = genesKeyInteractionPairSourceSink.get(gene);
		sinkInteraction = interactingGeneFrom.get(gene);
		if (sinkInteraction.contains(gene))
			sinkInteraction.remove(gene);
		List<String> sourceIneraction = new ArrayList<>();
		sourceIneraction = interactingGeneTo.get(gene);
		if (sourceIneraction.contains(gene))
			sourceIneraction.remove(gene);
		Boolean sink = false, source = false;
		String sinkGeneName = null, sourceGeneName = null;
		for (String sinkGene : sinkInteraction) {
			if (genesList.contains(sinkGene) && interctionTempTotally.contains(sinkGene)) {
				sink = true;
				sinkGeneName = sinkGene;
				break;
			}
		}
		for (String sourceGene : sourceIneraction) {
			if (genesList.contains(sourceGene) && interctionTempTotally.contains(sourceGene)) {
				source = true;
				sourceGeneName = sourceGene;
				break;
			}
		}
		if ((!(sink && source)) && (sinkGeneName != null || sourceGeneName != null)) {
			if (sink) {
				sinkCounter++;
				return "Sink";

			}
			if (source) {
				sourceCounter++;
				return "Source";
			}
		}
		return "Intermediate";

	}

	private static List<String> getUpstreamGenesNotinNetworkFile(String gene) {
		List<String> returnList = new ArrayList<>();
		List<String> interctionTempTotally = genesKeyInteractionPairSourceSink.get(gene);
		List<GenesInteractions> geneInteractionList = DataStore.getedgeClassListComplete();
		for (GenesInteractions edgeClassTemp : geneInteractionList) {
			if (edgeClassTemp.getColumn1().equalsIgnoreCase(gene)
					&& edgeClassTemp.getSymbol().equals(CustomEnum.reverseActivation)
					&& (!interctionTempTotally.contains(edgeClassTemp.getColmn2()))) {
				try {
					String column2Value = edgeClassTemp.getColmn2().toLowerCase();
					Double pvalueTemp = genesKeyPValuesPair.get(column2Value);
					if (!pvalueTemp.isNaN())
						if (pvalueTemp > 0 && pvalueTemp < 1) {
							returnList.add(edgeClassTemp.getColmn2());
						}
				} catch (Exception e) {
				}
			}
			if (edgeClassTemp.getColmn2().equalsIgnoreCase(gene)
					&& edgeClassTemp.getSymbol().equals(CustomEnum.activation)
					&& (!interctionTempTotally.contains(edgeClassTemp.getColumn1()))) {
				try {
					String column1Value = edgeClassTemp.getColumn1().toLowerCase();
					Double pvalueTemp = genesKeyPValuesPair.get(column1Value);
					if (!pvalueTemp.isNaN())
						if (pvalueTemp > 0 && pvalueTemp < 1) {
							returnList.add(edgeClassTemp.getColumn1().toLowerCase());
						}
				} catch (Exception e) {
				}
			}
		}
		return returnList;
	}

	private static List<String> getDownstreamGenesOnlyinNetworkFile(String gene) {
		List<String> returnList = new ArrayList<>();
		List<String> interctionTempTotally = genesKeyInteractionPairSourceSink.get(gene);
		List<GenesInteractions> geneInteractionList = DataStore.getedgeClassListComplete();
		for (GenesInteractions edgeClassTemp : geneInteractionList) {
			if (edgeClassTemp.getColumn1().equalsIgnoreCase(gene)
					&& edgeClassTemp.getSymbol().equals(CustomEnum.activation)
					&& (interctionTempTotally.contains(edgeClassTemp.getColmn2()))) {
				try {
					String column2Value = edgeClassTemp.getColmn2().toLowerCase();
					Double pvalueTemp = genesKeyPValuesPair.get(column2Value);
					if (!pvalueTemp.isNaN())
						if (pvalueTemp > 0 && pvalueTemp < 1) {
							returnList.add(edgeClassTemp.getColmn2());
						}
				} catch (Exception e) {
				}
			}
			if (edgeClassTemp.getColmn2().equalsIgnoreCase(gene)
					&& edgeClassTemp.getSymbol().equals(CustomEnum.reverseActivation)
					&& (interctionTempTotally.contains(edgeClassTemp.getColumn1()))) {
				try {
					String column1Value = edgeClassTemp.getColumn1().toLowerCase();
					Double pvalueTemp = genesKeyPValuesPair.get(column1Value);
					if (!pvalueTemp.isNaN())
						if (pvalueTemp > 0 && pvalueTemp < 1) {
							returnList.add(edgeClassTemp.getColumn1().toLowerCase());
						}
				} catch (Exception e) {
				}
			}
		}
		return returnList;
	}

	private static synchronized void writeInFIle(List<String> output) {
		Charset charset = StandardCharsets.UTF_8;
		try {
			Files.write(Paths.get(fileoutPathText), (output.toString() + "\n").getBytes(charset),
					StandardOpenOption.CREATE, StandardOpenOption.APPEND);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static synchronized void readNetworkFile(String fileName) {
		List<String> list = new ArrayList<>();
		try (Stream<String> stream = Files.lines(Paths.get(fileName))) {
			list = stream.map(String::toUpperCase).collect(Collectors.toList());
		} catch (IOException e) {
			e.printStackTrace();
		}
		List<String> uniqueEnteriesTemp = new ArrayList<>();
		List<GenesInteractions> genesInteractionsListCompleteTemp = new ArrayList<GenesInteractions>();
		for (String stringTemp : list) {
			if (stringTemp.equalsIgnoreCase(""))
				continue;
			String[] splitArray = stringTemp.split("\\s+");
			String from = splitArray[0].toLowerCase();
			String to = splitArray[1].toLowerCase();
			String symbol = splitArray[2].toLowerCase();
			if (!uniqueEnteriesTemp.contains(to))
				uniqueEnteriesTemp.add(to);
			if (!uniqueEnteriesTemp.contains(from))
				uniqueEnteriesTemp.add(from);
			GenesInteractions edgeClass = new GenesInteractions();
			edgeClass.setColumn1(from);
			edgeClass.setColmn2(to);
			edgeClass.setSymbol(symbol);
			genesInteractionsListCompleteTemp.add(edgeClass);
		}
		for (String uniqueKey : uniqueEnteriesTemp) {
			List<String> valuesList = new ArrayList<String>();
			for (GenesInteractions edgeClass : genesInteractionsListCompleteTemp) {
				if (edgeClass.getColumn1().equalsIgnoreCase(uniqueKey)) {
					valuesList.add(edgeClass.getColmn2());
				} else if (edgeClass.getColmn2().equalsIgnoreCase(uniqueKey)) {
					valuesList.add(edgeClass.getColumn1());
				}
			}
			genesKeyInteractionPairSourceSink.put(uniqueKey, valuesList);
		}
	}

	private static void refreshList() {
		List<GenesInteractions> genesInteractionsList = new ArrayList<>();
		for (GenesInteractions edgeClass : genesInteractionsListComplete) {
			String symbol = edgeClass.getSymbol();
			if (symbol.equalsIgnoreCase(CustomEnum.activation) || symbol.equalsIgnoreCase(CustomEnum.reverseActivation)
					|| symbol.equalsIgnoreCase(CustomEnum.Inhibition)
					|| symbol.equalsIgnoreCase(CustomEnum.reverseInhibitor)) {

				if (symbol.equalsIgnoreCase(CustomEnum.activation) || symbol.equalsIgnoreCase(CustomEnum.Inhibition))
					genesInteractionsList.add(edgeClass);
				else if (symbol.equalsIgnoreCase(CustomEnum.reverseActivation)
						|| symbol.equalsIgnoreCase(CustomEnum.reverseInhibitor)) {
					String firstColumn = edgeClass.getColumn1();
					edgeClass.setColumn1(edgeClass.getColmn2());
					edgeClass.setColmn2(firstColumn);
					genesInteractionsList.add(edgeClass);
				}
			}
		}
		genesInteractionsListComplete = new ArrayList<>();
		genesInteractionsListComplete = genesInteractionsList;
	}

	private static void createHashListFromEdegeClassList(int column) {
		refreshList();
		genesKeyInteractionPair = new HashMap<String, List<String>>();
		interactingGeneFrom = new HashMap<String, List<String>>();
		interactingGeneTo = new HashMap<String, List<String>>();
		int column1 = 0;
		int column2 = 0;
		for (String uniqueKey : uniqueGenesInInteraction) {
			List<String> valuesList = new ArrayList<String>();
			List<String> valuesList1 = new ArrayList<String>();
			List<String> valuesList2 = new ArrayList<String>();
			if (column == 0 || column == 2) {
				column1++;
				for (GenesInteractions edgeClass : genesInteractionsListComplete) {
					if (valuesList.size() > 1000)
						break;
					if (edgeClass.getColumn1().equalsIgnoreCase(uniqueKey)) {
						try {
							String column2Value = edgeClass.getColmn2().toLowerCase();
							Double pvalueTemp = genesKeyPValuesPair.get(column2Value);
							if (!pvalueTemp.isNaN())
								if (pvalueTemp > 0 && pvalueTemp < 1) {
									valuesList.add(column2Value);
									valuesList1.add(column2Value);
								}
						} catch (Exception e) {
						}
					}
				}
				interactingGeneTo.put(uniqueKey.toLowerCase(), valuesList1);
			}
			if (column == 0 || column == 1) {
				column2++;
				for (GenesInteractions edgeClass : genesInteractionsListComplete) {
					if (valuesList.size() > 1000)
						break;
					if (edgeClass.getColmn2().equalsIgnoreCase(uniqueKey)) {
						try {
							String column1Value = edgeClass.getColumn1().toLowerCase();
							Double pvalueTemp = genesKeyPValuesPair.get(column1Value);
							if (!pvalueTemp.isNaN())
								if (pvalueTemp > 0 && pvalueTemp < 1) {
									valuesList.add(column1Value);
									valuesList2.add(column1Value);
								}
						} catch (Exception e) {
						}
					}
				}
				interactingGeneFrom.put(uniqueKey.toLowerCase(), valuesList2);
			}
			genesKeyInteractionPair.put(uniqueKey.toLowerCase(), valuesList);
		}
	}
}