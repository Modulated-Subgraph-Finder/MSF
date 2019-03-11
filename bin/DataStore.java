import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.io.BufferedWriter;
import java.io.FileWriter;

import javax.print.attribute.standard.MediaSize.NA;
import javax.print.attribute.standard.MediaSize.Other;


public class DataStore {

	public static GenesInfo geneInfo = new GenesInfo();
	public static List<GenesInfo> geneInfoList = new ArrayList<GenesInfo>();
	public static Map<String, Double> genesKeyPValuesPair = new HashMap<String, Double>();
	public static Map<String, String> genesKeyLogFoldChangesPair = new HashMap<String, String>();	
	public static List<String> uniqueGenesInInteraction = new ArrayList<String>();
	public static List<GenesInteractions> genesInteractionsListComplete = new ArrayList<GenesInteractions>();
	public static List<GenesInteractions> genesInteractionsListCompleteTemp = new ArrayList<GenesInteractions>();
	static Map<String, List<TempModel>> interactingGeneFrom = new HashMap<String, List<TempModel>>();
	static Map<String, List<TempModel>> interactingGeneTo = new HashMap<String, List<TempModel>>();
	static Map<String, List<String>> genesKeyInteractionPair = new HashMap<String, List<String>>();

	static Map<String, List<TempModel>> genesKeyInteractionPairTemp = new HashMap<String, List<TempModel>>();

	static Hashtable<String, String> geneInfoHashTable = new Hashtable<>();
	

	static String outputPath;
	static String dEgType;
	private static DataStore instance = null;

	public static String getOutputPath() {
		return outputPath;
	}

	public static void setOutputPath(String outputPath) {
		DataStore.outputPath = outputPath;
	}

	public static void DataStor(String psheet, String interaction, String outputPathX, String dEgtype) {
		outputPath = outputPathX;
		dEgType = dEgtype;
		readPValuesFromTextFile(psheet);
		readInteractionsFromTextFile(interaction);
		createHashListFromEdegeClassList();
		System.out.println("Finding Initial Graphs");
	}

	public static void main(String[] args) throws Exception {
		{

		}
	}

	public static DataStore getInstance() {
		if (instance == null) {
			instance = new DataStore();
		}
		return instance;
	}

	protected static GenesInfo getpSheet() {
		return geneInfo;
	}

	protected static List<GenesInfo> getpSheetList() {
		return geneInfoList;
	}

	protected static Map<String, Double> getpsheetKeyValyePair() {
		return genesKeyPValuesPair;
	}

	protected static Map<String, String> getpsheetKeyFolderValue() {
		return genesKeyLogFoldChangesPair;
	}
	protected static List<String> getedgeListSet() {
		return uniqueGenesInInteraction;
	}

	protected static List<GenesInteractions> getedgeClassListComplete() {
		return genesInteractionsListComplete;
	}

	protected static Map<String, List<String>> getedgeListHashMap() {
		return genesKeyInteractionPair;
	}

	protected static Hashtable<String, String> getgeneInfoHashTable() {
		return geneInfoHashTable;
	}

	public static void readPValuesFromTextFile(String fileName) {

		List<String> list = new ArrayList<>();
		try (Stream<String> stream = Files.lines(Paths.get(fileName))) {
			list = stream.map(String::toUpperCase).collect(Collectors.toList());
		} catch (IOException e) {
			e.printStackTrace();
		}
		int count = 0;
		list.remove(0);
		for (String string : list) { 
			count = count + 1;
			geneInfo = new GenesInfo();
			String[] splitArray = string.split("\\s+"); 
			String geneName = splitArray[0].replace("\"", ""); 
			String pvalue = "";
			String logFoldChange = "";
		
			if (dEgType.equals("DEseq2")){
				try{
				pvalue = splitArray[5];
				logFoldChange = splitArray[2];

				} catch (Exception e) {
					System.out.println("Please give Tab-seperated file");
				}
			}
			else {
				try{
					pvalue = splitArray[4];
					logFoldChange = splitArray[1];
					} catch (Exception e) {
						System.out.println("Please give Tab-seperated file");
					}
			}
		
				geneInfo.setGeneName(geneName.toLowerCase());
			if (pvalue.equalsIgnoreCase("NA"))
				continue;
			double tempasd = (Double.parseDouble(pvalue));
			if (tempasd == 2D)
				tempasd = 1;
			geneInfo.setpValue(tempasd);
			if (geneInfo.getpValue() > 0 && geneInfo.getpValue() < 1)
				geneInfoList.add(geneInfo);
			genesKeyPValuesPair.put(geneInfo.getGeneName().toLowerCase(), geneInfo.getpValue());
			genesKeyLogFoldChangesPair.put(geneInfo.getGeneName().toLowerCase(), logFoldChange);
			geneInfoHashTable.put(geneName.toLowerCase(), String.valueOf(splitArray[1]));
		}
	}

	public static void readInteractionsFromTextFile(String fileName) {
		List<String> list = new ArrayList<>();
		try (Stream<String> stream = Files.lines(Paths.get(fileName))) {
			list = stream.map(String::toUpperCase).collect(Collectors.toList());
		} catch (IOException e) {
			e.printStackTrace();
		}
		int count = 0;
		int mismatchedCount = 0;
		int listSize = list.size();
		for (String string : list) {
			count = count + 1;
			GenesInteractions edgeClass = new GenesInteractions();
			String[] splitArray = string.split("\\s+");
			String from = splitArray[0].toLowerCase();
			String to = splitArray[1].toLowerCase();
			String symbol = splitArray[2];
			Boolean pvalueFound = true;
			try {
				Double pvalue = genesKeyPValuesPair.get(from);
				Double pvalue2 = genesKeyPValuesPair.get(to);
				if (pvalue == null || pvalue == Double.NaN || pvalue2 == null || pvalue2 == Double.NaN) 
																									
				{
					mismatchedCount ++;
					continue;
				}
			} catch (Exception e) {
				pvalueFound = false;
			}
			if (!pvalueFound)
				continue;
			if (!uniqueGenesInInteraction.contains(from))
				uniqueGenesInInteraction.add(from);
			if (!uniqueGenesInInteraction.contains(to))
				uniqueGenesInInteraction.add(to);
			edgeClass.setColumn1(from);
			edgeClass.setColmn2(to);
			edgeClass.setSymbol(symbol);
			genesInteractionsListComplete.add(edgeClass);
		}
		{
		if (mismatchedCount==listSize){
			System.out.println("Gene Identifers do not match in the input files.");
		System.exit(0);}
		}
	}

	static void createHashListFromEdegeClassList() {
		for (String uniqueKey : uniqueGenesInInteraction) {
			List<String> valuesList = new ArrayList<String>();
			for (GenesInteractions edgeClass : genesInteractionsListComplete) {
				if (edgeClass.getColumn1().equalsIgnoreCase(uniqueKey)) {
					valuesList.add(edgeClass.getColmn2());
				}
			}
			for (GenesInteractions edgeClass : genesInteractionsListComplete) {
				if (edgeClass.getColmn2().equalsIgnoreCase(uniqueKey)) {
					valuesList.add(edgeClass.getColumn1());
				}
			}
			genesKeyInteractionPair.put(uniqueKey, valuesList);
		}
	}
	

	private static void refreshList() {
		List<GenesInteractions> genesInteractionsList = new ArrayList<>();
		for (GenesInteractions edgeClass : genesInteractionsListComplete) {
			String symbol = edgeClass.getSymbol();
				genesInteractionsList.add(edgeClass);

		}
		genesInteractionsListCompleteTemp = new ArrayList<>();
		genesInteractionsListCompleteTemp = genesInteractionsList;
	}

	public static Map<String, List<TempModel>> getInteractingGeneFrom() {
		return interactingGeneFrom;
	}

	public static void setInteractingGeneFrom(Map<String, List<TempModel>> interactingGeneFrom) {
		DataStore.interactingGeneFrom = interactingGeneFrom;
	}

	public static Map<String, List<TempModel>> getInteractingGeneTo() {
		return interactingGeneTo;
	}

	public static void setInteractingGeneTo(Map<String, List<TempModel>> interactingGeneTo) {
		DataStore.interactingGeneTo = interactingGeneTo;
	}

	static void createHashListFromEdegeClassListTemp(int column) {
		refreshList();
		interactingGeneFrom = new HashMap<String, List<TempModel>>();
		interactingGeneTo = new HashMap<String, List<TempModel>>();
		int column1 = 0;
		int column2 = 0;
		for (String uniqueKey : uniqueGenesInInteraction) {
			List<TempModel> valuesList = new ArrayList<TempModel>();
			List<TempModel> valuesList1 = new ArrayList<TempModel>();
			List<TempModel> valuesList2 = new ArrayList<TempModel>();
			if (column == 0 || column == 2) {
				column1++;
				for (GenesInteractions edgeClass : genesInteractionsListCompleteTemp) {
					if (valuesList.size() > 1000)
						break;
					if (edgeClass.getColumn1().equalsIgnoreCase(uniqueKey)) {
						try {
							String column2Value = edgeClass.getColmn2().toLowerCase();
							Double pvalueTemp = genesKeyPValuesPair.get(column2Value);
							if (!pvalueTemp.isNaN())
								if (pvalueTemp > 0 && pvalueTemp < 1) {
									TempModel tempModel = new TempModel();
									tempModel.setInteractingGene(column2Value.trim());
									tempModel.setSymbol(edgeClass.getSymbol());
									valuesList.add(tempModel);
									valuesList1.add(tempModel);
								}
						} catch (Exception e) {
						}
					}
				}
				interactingGeneTo.put(uniqueKey.toLowerCase().trim(), valuesList1);
			} 
			if (column == 0 || column == 1) {
				column2++;
				for (GenesInteractions edgeClass : genesInteractionsListCompleteTemp) {
					if (valuesList.size() > 1000)
						break;
					if (edgeClass.getColmn2().equalsIgnoreCase(uniqueKey)) {
						try {
						String column1Value = edgeClass.getColumn1().toLowerCase();
							Double pvalueTemp = genesKeyPValuesPair.get(column1Value);
							if (!pvalueTemp.isNaN())
								if (pvalueTemp > 0 && pvalueTemp < 1) {
									TempModel tempModel = new TempModel();
									tempModel.setInteractingGene(column1Value.trim());
									tempModel.setSymbol(edgeClass.getSymbol());
									valuesList.add(tempModel);
									valuesList2.add(tempModel);
								}
						} catch (Exception e) {
						}
					}
				}
				interactingGeneFrom.put(uniqueKey.toLowerCase().trim(), valuesList2);
			}
			genesKeyInteractionPairTemp.put(uniqueKey.toLowerCase().trim(), valuesList);
		}
	}
	
}