import java.io.BufferedWriter;
import java.io.FileInputStream;
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
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class ExtensionMerging {

	static Integer interactingListCounter = 0;
	static Integer interactingListCounterCheck = 0;
	static List<Integer> alreadyIteratedPath = new ArrayList<>();
	static List<String> finalizingPath = new ArrayList<>();
	static List<String> step2GnesDone = new ArrayList<>();
	static List<String> genesDone = new ArrayList<>();
	static Set<String> aaa = new HashSet<String>();
	public static List<ArrayList<String>> AllMergedpaths = new ArrayList<ArrayList<String>>();
	static Map<String, List<TempModel>> interactingGeneFrom = new HashMap<String, List<TempModel>>();
	static Map<String, List<TempModel>> interactingGeneTo = new HashMap<String, List<TempModel>>();
	static String Step2Output;
	static String Step3Output;
	static List<GenesInfo> geneInfoList = new ArrayList<GenesInfo>();
	static Map<String, Double> genesKeyPValuesPair = new HashMap<String, Double>();
	static List<GenesInteractions> genesInteractionsListComplete = new ArrayList<GenesInteractions>();
	static List<String> uniqueGenesInInteraction = new ArrayList<String>();
	static GenesInfo genesInfo = new GenesInfo();
	static Map<String, List<String>> genesKeyInteractionPair = new HashMap<String, List<String>>();
	static long THRESHOLD = 2L;
	static int rowCounter = 0;
	static List<MergingPojo> mergedPath = new ArrayList<>();
	static List<Integer> mergedPathPathNumber = new ArrayList<>();
	static Hashtable<Integer, List<String>> hashTableList = new Hashtable<Integer, List<String>>();
	static FileInputStream fileInputSteam;
	static List<StringListSorting> sortingList;
	static String network_File, mergingNetwork;
	static int arrayListSize = 0, extensionnLimit = 2, mergingLimit = 1;
	static Boolean printExtraFiles = false;

	public static void main(List<ArrayList<String>> Allpaths, String filePath, int eLimit, int mLimit,
			Boolean printExtraFilesTemp) throws Exception {
		extensionnLimit = eLimit;
		mergingLimit = mLimit;
		System.out.println("Extending and Merging Graphs");
		Step2Output = DataStore.getOutputPath() + "ExtendedGraphs.text";
		Step3Output = DataStore.getOutputPath() + "MergedGraphs.text";
		network_File = filePath + "NetworkFile.text";
		mergingNetwork = filePath + "merging.text";
		printExtraFiles = printExtraFilesTemp;
		genesInfo = DataStore.getpSheet();
		geneInfoList = DataStore.getpSheetList();
		genesInteractionsListComplete = DataStore.getedgeClassListComplete();
		uniqueGenesInInteraction = DataStore.getedgeListSet();
		genesKeyPValuesPair = DataStore.getpsheetKeyValyePair();
		genesKeyInteractionPair = DataStore.getedgeListHashMap();
		arrayListSize = Allpaths.size();
		int listNumer = 1;
		for (ArrayList<String> arrayList : Allpaths) {
			setGenesList(arrayList, (listNumer));
			listNumer++;
			for (String string : arrayList) {
				step2GnesDone.add(string);
			}
		}
		extendingPaths();
		try {
			Thread.sleep(3000);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		List<String> tempStringList = null;
		for (int i = 1; i <= arrayListSize; i++) {
			tempStringList = getGenesList(i);
			Integer currentNumber = i;
			if (alreadyIteratedPath.contains(currentNumber))
				continue;
			List<String> sending = tempStringList;
			try {
				interactingListCounter = 0;
				interactingListCounterCheck = 0;
				interactingListCounterCheck = sending.size();
				for (interactingListCounter = 0; interactingListCounter < interactingListCounterCheck; interactingListCounter++) {
					String gene = "";
					gene = sending.get(interactingListCounter);
					mergingPaths(gene, i);
				}
			} catch (Exception e) {
			}
		}
		for (int x = 1; x <= arrayListSize; x++) {
			if (getGenesList(x).size() > 1)
				AllMergedpaths.add((ArrayList<String>) getGenesList(x));
		}

		if (printExtraFiles == true) {
			writeInFIleMergedPaths();
		}
		System.out.println("Extension and Merging Done");
		FinalThreads.main(AllMergedpaths,network_File);
	}

	private static synchronized void writeInFIleExtension() {

		Charset charset = StandardCharsets.UTF_8;
		try {
			for (int xc = 1; xc <= arrayListSize; xc++) {
				List<String> tempOrg = getGenesList(xc);
				List<Double> pavlueListString = new ArrayList<>();
				pavlueListString.clear();
				for (String string : tempOrg) {
					pavlueListString.add(genesKeyPValuesPair.get(string));
				}
				Files.write(Paths.get(Step2Output),
						(tempOrg.toString() + " " + GenericFunctions.hartungFunction(pavlueListString) + "\n")
								.getBytes(charset),
						StandardOpenOption.CREATE, StandardOpenOption.APPEND);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static void extendingPaths() {
		Set<Integer> pathSet = new HashSet<Integer>();
		interactingGeneFrom = DataStore.getInteractingGeneFrom();
		interactingGeneTo = DataStore.getInteractingGeneTo();
		HashMap<String, String> geneInteractingAncestor = new HashMap<>();
		for (int pathNumber = 1; pathNumber <= arrayListSize; pathNumber++) {
			List<String> tempList = getGenesList(pathNumber);
			Boolean pathExtended = false;
			for (String string : tempList) {
				List<String> iteratingGenes = new ArrayList<>();
				while (true) {
					List<String> threeLengthPath = genesPathSize3(string, iteratingGenes);
					if (threeLengthPath == null || threeLengthPath.size() < 1)
						break;
					List<Double> extendedPathCP = new ArrayList<>();
					List<Double> originalPathCP = new ArrayList<>();
					List<String> newExtendedPath = new ArrayList<>();
					List<String> tempOrg = getGenesList(pathNumber);
					for (String gene : tempOrg) {
						newExtendedPath.add(gene);
						originalPathCP.add(genesKeyPValuesPair.get(gene));
						extendedPathCP.add(genesKeyPValuesPair.get(gene));
					}
					for (String string1 : threeLengthPath) {
						newExtendedPath.add(string1);
						extendedPathCP.add(genesKeyPValuesPair.get(string1));
					}
					Double originalPathHartung = GenericFunctions.hartungFunction(originalPathCP);
					if (originalPathHartung == 0.0)
						originalPathHartung = 2E-16;
					Double extendedPathHartung = GenericFunctions.hartungFunction(extendedPathCP);
					if (originalPathHartung > extendedPathHartung) {
						for (int countX = 0; countX < threeLengthPath.size(); countX++) {
							String ax = threeLengthPath.get(countX);
							if (countX == 0)
								geneInteractingAncestor.put(ax, string.trim());
							else {
								geneInteractingAncestor.put(threeLengthPath.get(countX),
										threeLengthPath.get(countX - 1));
							}
						}
						for (String as : threeLengthPath) {
							step2GnesDone.add(as);
							String parentGeneOFChildGene = geneInteractingAncestor.get(as);
							List<TempModel> tempModelList = interactingGeneTo.get(parentGeneOFChildGene);
							List<TempModel> xx = tempModelList.stream()
									.filter(x -> x.interactingGene.equalsIgnoreCase(as.trim()))
									.collect(Collectors.toList());
							TempModel tempp = new TempModel();
							ArrayList<String> printa = new ArrayList<>();
							String printx = "";
							if (xx.size() < 1) {
								tempModelList = interactingGeneFrom.get(parentGeneOFChildGene);
								xx = tempModelList.stream().filter(x -> x.interactingGene.equalsIgnoreCase(as.trim()))
										.collect(Collectors.toList());
								if (xx.size() > 0) {
									tempp = xx.get(0);
									printx = tempp.getInteractingGene() + " " + parentGeneOFChildGene + " "
											+ tempp.getSymbol();
									if (tempp.getSymbol().equals(CustomEnum.activation))
										printx = printx + " 1 "+ " 0 ";
									else if (tempp.getSymbol().equals(CustomEnum.Inhibition))
										printx = printx + " 1 "+ " 0 ";
									else if (tempp.getSymbol().equals(CustomEnum.reverseInhibitor))
										printx = printx + " 0 "+ " 1 ";
									else if (tempp.getSymbol().equals(CustomEnum.reverseActivation))
										printx = printx + " 0 "+ " 1 ";
									else if (tempp.getSymbol().equals(CustomEnum.inhibitionActivation)|| tempp.getSymbol().equals(CustomEnum.inhibitionInhibition) || tempp.getSymbol().equals(CustomEnum.activationInhibition)|| tempp.getSymbol().equals(CustomEnum.activationActivation))
										printx = printx + " 1 "+ " 1 ";
								}
							} else {
								tempp = xx.get(0);
								printx = parentGeneOFChildGene + " " + tempp.getInteractingGene() + " "
										+ tempp.getSymbol();
								if (tempp.getSymbol().equals(CustomEnum.activation))
									printx = printx + " 1 "+ " 0 ";
								else if (tempp.getSymbol().equals(CustomEnum.Inhibition))
									printx = printx + " 1 "+ " 0 ";
								else if (tempp.getSymbol().equals(CustomEnum.reverseInhibitor))
									printx = printx + " 0 "+ " 1 ";
								else if (tempp.getSymbol().equals(CustomEnum.reverseActivation))
									printx = printx + " 0 "+ " 1 ";
								else if (tempp.getSymbol().equals(CustomEnum.inhibitionActivation)|| tempp.getSymbol().equals(CustomEnum.inhibitionInhibition) || tempp.getSymbol().equals(CustomEnum.activationInhibition)|| tempp.getSymbol().equals(CustomEnum.activationActivation))
									printx = printx + " 1 "+ " 1 ";
							}
							writeInFIle_NetworkFile(printx);
							// close file writing
						}
						pathSet.add(pathNumber);
						setGenesList(newExtendedPath, pathNumber);
						pathExtended = true;
					} else if (checkPathwithSig(threeLengthPath, pathNumber, string)) {
						pathExtended = true;
					}
				}
			}

		}
		if (printExtraFiles == true) {
			writeInFIleExtension();
		}
	}

	private static synchronized void writeInFIle_NetworkFile(String output) {

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

	private static boolean checkPathwithSig(List<String> ext, int pathNumber, String string) {
		boolean upadteCheck = false;
		for (int qw = ext.size(); qw >= 1; qw--) {
			List<String> temp = new ArrayList<String>();
			for (int xc = 0; xc < qw - 1; xc++) {
				temp.add(ext.get(xc));
			}
			List<Double> originalPathCP = new ArrayList<>();
			List<Double> extendedPathCP = new ArrayList<>();
			List<String> newExtendedPath = new ArrayList<>();
			List<String> tempOrg = getGenesList(pathNumber);
			for (String gene : tempOrg) {
				newExtendedPath.add(gene);
				originalPathCP.add(genesKeyPValuesPair.get(gene));
				extendedPathCP.add(genesKeyPValuesPair.get(gene));
			}
			for (String string1 : temp) {
				newExtendedPath.add(string1);
				extendedPathCP.add(genesKeyPValuesPair.get(string1));
			}
			Double originalPathHartung = GenericFunctions.hartungFunction(originalPathCP);
			if (originalPathHartung == 0.0)
				originalPathHartung = 2E-16;
			Double extendedPathHartung = GenericFunctions.hartungFunction(extendedPathCP);
			if (originalPathHartung > extendedPathHartung) {
				for (String as : temp) {
					step2GnesDone.add(as);
				}
				HashMap<String, String> geneInteractingAncestor = new HashMap<>();
				for (int countX = 0; countX < temp.size(); countX++) {
					String ax = temp.get(countX);
					if (countX == 0)
						geneInteractingAncestor.put(ax, string.trim());
					else {
						geneInteractingAncestor.put(temp.get(countX), temp.get(countX - 1));
					}
				}

				for (String as : temp) {
					step2GnesDone.add(as);
					String parentGeneOFChildGene = geneInteractingAncestor.get(as);
					List<TempModel> tempModelList = interactingGeneTo.get(parentGeneOFChildGene);
					List<TempModel> xx = tempModelList.stream()
							.filter(x -> x.interactingGene.equalsIgnoreCase(as.trim())).collect(Collectors.toList());
					TempModel tempp = new TempModel();
					ArrayList<String> printa = new ArrayList<>();
					String printx = "";
					if (xx.size() < 1) {
						tempModelList = interactingGeneFrom.get(parentGeneOFChildGene);
						xx = tempModelList.stream().filter(x -> x.interactingGene.equalsIgnoreCase(as.trim()))
								.collect(Collectors.toList());
						if (xx.size() > 0) {
							tempp = xx.get(0);
							printx = tempp.getInteractingGene() + " " + parentGeneOFChildGene + " " + tempp.getSymbol();
							if (tempp.getSymbol().equals(CustomEnum.activation))
								printx = printx + " 1 "+ " 0 ";
							else if (tempp.getSymbol().equals(CustomEnum.Inhibition))
								printx = printx + " 1 "+ " 0 ";
							else if (tempp.getSymbol().equals(CustomEnum.reverseInhibitor))
								printx = printx + " 0 "+ " 1 ";
							else if (tempp.getSymbol().equals(CustomEnum.reverseActivation))
								printx = printx + " 0 "+ " 1 ";
							else if (tempp.getSymbol().equals(CustomEnum.inhibitionActivation)|| tempp.getSymbol().equals(CustomEnum.inhibitionInhibition) || tempp.getSymbol().equals(CustomEnum.activationInhibition)|| tempp.getSymbol().equals(CustomEnum.activationActivation))
								printx = printx + " 1 "+ " 1 ";
						}
					} else {
						tempp = xx.get(0);
						printx = parentGeneOFChildGene + " " + tempp.getInteractingGene() + " " + tempp.getSymbol();
						if (tempp.getSymbol().equals(CustomEnum.activation))
							printx = printx + " 1 "+ " 0 ";
						else if (tempp.getSymbol().equals(CustomEnum.Inhibition))
							printx = printx + " 1 "+ " 0 ";
						else if (tempp.getSymbol().equals(CustomEnum.reverseInhibitor))
							printx = printx + " 0 "+ " 1 ";
						else if (tempp.getSymbol().equals(CustomEnum.reverseActivation))
							printx = printx + " 0 "+ " 1 ";
						else if (tempp.getSymbol().equals(CustomEnum.inhibitionActivation)|| tempp.getSymbol().equals(CustomEnum.inhibitionInhibition) || tempp.getSymbol().equals(CustomEnum.activationInhibition)|| tempp.getSymbol().equals(CustomEnum.activationActivation))
							printx = printx + " 1 "+ " 1 ";
					}
					writeInFIle_NetworkFile(printx);
					// close file writing
				}

				// clsoing here
				setGenesList(newExtendedPath, pathNumber);
				upadteCheck = true;
				return true;

			}
		}
		return upadteCheck;

	}

	private static List<String> genesPathSize3(String gene, List<String> genesListIterating) {
		List<String> returnString = new ArrayList<String>();
		String int1 = smallestPvalueInteractingGene(gene, returnString, genesListIterating);
		if (int1 != null) {
			returnString.add(int1);
			if (extensionnLimit == 1)
				return returnString;
			String int2 = smallestPvalueInteractingGene(int1, returnString, genesListIterating);
			if (int2 != null) {
				returnString.add(int2);
				if (extensionnLimit == 2)
					return returnString;
				String int3 = smallestPvalueInteractingGene(int2, returnString, genesListIterating);
				if (int3 != null) {
					returnString.add(int3);

					if (extensionnLimit == 3)
						return returnString;
				}
			}
		}
		return returnString;
	}

	private static String smallestPvalueInteractingGene(String gene, List<String> temp, List<String> itertingGenes) {
		List<String> interction = genesKeyInteractionPair.get(gene);
		interction = sortGeneInfoList(interction);
		for (String intr1 : interction) {
			if ((!step2GnesDone.contains(intr1)) && (!temp.contains(intr1)) && (!itertingGenes.contains(intr1))) {
				try {
					itertingGenes.add(intr1);
					Double pvalue = genesKeyPValuesPair.get(intr1);
					if (pvalue < 1)
						return intr1;
				} catch (Exception e) {

				}
			}
		}
		return null;
	}

	private static void mergingPaths(String gene, int listNumber) throws InterruptedException {
		HashMap<String, String> geneInteractingAncestor = new HashMap<>();
		Set<String> iterationSet = new HashSet<String>();
		iterationSet.add(gene);
		double keyPValue1 = genesKeyPValuesPair.get(gene);
		List<String> interactingGenes1 = genesKeyInteractionPair.get(gene);
		interactingGenes1 = sortGeneInfoList(interactingGenes1);
		List<String> parentGene = new ArrayList<>();
		parentGene.add(gene);
		geneInteractingAncestor.put(gene, "-");
		checkIfIneractionExist(gene, interactingGenes1, listNumber, 1, parentGene, geneInteractingAncestor);
		if (interactingGenes1 != null)
			for (String item1 : interactingGenes1) {
				if (!step2GnesDone.contains(item1))
					if (!iterationSet.contains(item1)) {
						iterationSet.add(item1);
						parentGene = new ArrayList<>();
						List<String> interactingGenes2 = genesKeyInteractionPair.get(item1);
						interactingGenes2 = sortGeneInfoList(interactingGenes2);
						parentGene.add(gene);
						geneInteractingAncestor.put(item1, gene);
						parentGene.add(item1);
						checkIfIneractionExist(item1, interactingGenes2, listNumber, 2, parentGene,geneInteractingAncestor);
						for (String item2 : interactingGenes2) {
							if (!step2GnesDone.contains(item2))
								if (!iterationSet.contains(item2)) {
									iterationSet.add(item2);
									parentGene = new ArrayList<>();
									List<String> interactingGenes3 = genesKeyInteractionPair.get(item2);
									interactingGenes3 = sortGeneInfoList(interactingGenes3);
									parentGene.add(gene);
									parentGene.add(item1);
									geneInteractingAncestor.put(item2, item1);
									parentGene.add(item2);
									if (mergingLimit > 1)
										// return;
										{
									checkIfIneractionExist(item2, interactingGenes3, listNumber, 3, parentGene,
											geneInteractingAncestor);
									for (String item3 : interactingGenes3) {
										if (!step2GnesDone.contains(item3))
											if (!iterationSet.contains(item3)) {
												iterationSet.add(item3);
												parentGene = new ArrayList<>();
												List<String> interactingGenes4 = genesKeyInteractionPair.get(item3);
												interactingGenes4 = sortGeneInfoList(interactingGenes4);
												parentGene.add(gene);
												parentGene.add(item1);
												parentGene.add(item2);
												geneInteractingAncestor.put(item3, item2);
												parentGene.add(item3);
												if (mergingLimit > 2)
												
												checkIfIneractionExist(item3, interactingGenes4, listNumber, 4,
														parentGene, geneInteractingAncestor);
											}
									}
								}// first new bracket
								}
						}
					}
			}
	}

	private static void checkIfIneractionExist(String genedd, List<String> interactingGenes1, int listNumber, int level,
			List<String> parentsGenes, HashMap<String, String> geneInteractingAncestor) throws InterruptedException {
		if (parentsGenes.size() < 2)
			return;
		if (interactingGenes1 != null)
			for (String string : interactingGenes1) {
				List<String> inerationWith = null;
				for (int x = 1; x <= arrayListSize; x++) {
					if (x == listNumber)
						continue;
					inerationWith = getGenesList(x);
					if (inerationWith == null)
						continue;
					Integer asd = x;
					if (inerationWith.contains(string)) {
						boolean geneExist = false;
						List<String> checkList = new ArrayList<String>(parentsGenes);
						checkList.remove(0);
						for (String string2 : checkList) {
							if (inerationWith.contains(string2)) {
								geneExist = true;
							}
						}
						if (geneExist)
							return;
						List<String> originlParentList = new ArrayList<>();
						originlParentList = getGenesList(listNumber);
						if (originlParentList == null)
							continue;
						boolean geneExist1 = false;
						List<String> checkList1 = new ArrayList<String>(parentsGenes);
						checkList1.remove(0);
						for (String string2 : checkList1) {
							if (originlParentList.contains(string2)) {
								geneExist1 = true;
							}
						}
						if (geneExist1)
							return;
						List<Double> originalParentPValue = new ArrayList<>();
						for (String pGene : originlParentList) {
							originalParentPValue.add(genesKeyPValuesPair.get(pGene));
						}
						Double OriginlHartungResults = GenericFunctions.hartungFunction(originalParentPValue);
						List<Double> originalInteractionPValue = new ArrayList<>();
						for (String pGene : inerationWith) {
							originalInteractionPValue.add(genesKeyPValuesPair.get(pGene));
						}
						Double interactionHartungResults = GenericFunctions.hartungFunction(originalInteractionPValue);
						List<String> newExtendedPath = new ArrayList<>();
						List<Double> newExtendedPathPValue = new ArrayList<>();
						Double hartungResultNew = 0.0;
						for (String s : originlParentList) {
							newExtendedPath.add(s);
							newExtendedPathPValue.add(genesKeyPValuesPair.get(s));
						}
						for (String s : checkList) {
							newExtendedPath.add(s);
							newExtendedPathPValue.add(genesKeyPValuesPair.get(s));
						}
						for (String s : inerationWith) {
							newExtendedPath.add(s);
							newExtendedPathPValue.add(genesKeyPValuesPair.get(s));
						}

						hartungResultNew = GenericFunctions.hartungFunction(newExtendedPathPValue);
						Boolean resultLessParent = (OriginlHartungResults > hartungResultNew);
						Boolean resultLessExtended = (interactionHartungResults > hartungResultNew);
						if (resultLessParent && resultLessExtended) {
							for (String sqw : checkList) {
								step2GnesDone.add(sqw);
								printingFunction(sqw, geneInteractingAncestor);
							}
							geneInteractingAncestor.put(string, genedd);
							printingFunction(string, geneInteractingAncestor);
							aaa.add(Integer.toString(listNumber));
							aaa.add(Integer.toString(x));
							MergingPojo merged = new MergingPojo();
							merged.setMergedPath(newExtendedPath);
							merged.setPathNumber_1(listNumber);
							setGenesList(newExtendedPath, listNumber);
							setGenesList(new ArrayList<String>(), x);
							interactingListCounter = 0;
							interactingListCounterCheck = newExtendedPath.size();
							merged.setPathNumber_2(x);
							merged.setCombineCP(hartungResultNew);
							alreadyIteratedPath.add(listNumber);
							alreadyIteratedPath.add(x);
							for (String doubleString : newExtendedPath) {
								if (!finalizingPath.contains(doubleString))
									finalizingPath.add(doubleString);
							}
							Boolean found = false;
							for (int check = 0; check < mergedPath.size(); check++) {
								MergingPojo merged1 = mergedPath.get(check);
								Boolean matchedList = merged.compareMergingPojosObject(merged1);
								if (matchedList) {
									found = true;
									if (hartungResultNew < merged1.getCombineCP()) {
										mergedPath.remove(check);
										mergedPath.add(check, merged);
									}
									break;
								}

							}
							if (!found) {
								mergedPathPathNumber.add(listNumber);
								mergedPathPathNumber.add(x);
								mergedPath.add(merged);
							}
						}

					}
				}
			}
	}

	private static void printingFunction(String gene, HashMap<String, String> geneInteractingAncestor) {
		String parentGeneOFChildGene = geneInteractingAncestor.get(gene);
		List<TempModel> tempModelList = interactingGeneTo.get(parentGeneOFChildGene);
		List<TempModel> xx = tempModelList.stream().filter(x1 -> x1.interactingGene.equalsIgnoreCase(gene))
				.collect(Collectors.toList());
		TempModel tempp = new TempModel();
		ArrayList<String> printa = new ArrayList<>();
		String printx = "";
		if (xx.size() < 1) {
			tempModelList = interactingGeneFrom.get(parentGeneOFChildGene);
			xx = tempModelList.stream().filter(x2 -> x2.interactingGene.equalsIgnoreCase(gene))
					.collect(Collectors.toList());
			if (xx.size() > 0) {
				tempp = xx.get(0);
				printx = tempp.getInteractingGene() + " " + parentGeneOFChildGene + " " + tempp.getSymbol();
				if (tempp.getSymbol().equals(CustomEnum.activation))
					printx = printx + " 1 "+ " 0 ";
				else if (tempp.getSymbol().equals(CustomEnum.Inhibition))
					printx = printx + " 1 "+ " 0 ";
				else if (tempp.getSymbol().equals(CustomEnum.reverseInhibitor))
					printx = printx + " 0 "+ " 1 ";
				else if (tempp.getSymbol().equals(CustomEnum.reverseActivation))
					printx = printx + " 0 "+ " 1 ";
				else if (tempp.getSymbol().equals(CustomEnum.inhibitionActivation)|| tempp.getSymbol().equals(CustomEnum.inhibitionInhibition) || tempp.getSymbol().equals(CustomEnum.activationInhibition)|| tempp.getSymbol().equals(CustomEnum.activationActivation))
					printx = printx + " 1 "+ " 1 ";
			}

		} else {
			tempp = xx.get(0);
			printx = parentGeneOFChildGene + " " + tempp.getInteractingGene() + " " + tempp.getSymbol();
			if (tempp.getSymbol().equals(CustomEnum.activation))
				printx = printx + " 1 "+ " 0 ";
			else if (tempp.getSymbol().equals(CustomEnum.Inhibition))
				printx = printx + " 1 "+ " 0 ";
			else if (tempp.getSymbol().equals(CustomEnum.reverseInhibitor))
				printx = printx + " 0 "+ " 1 ";
			else if (tempp.getSymbol().equals(CustomEnum.reverseActivation))
				printx = printx + " 0 "+ " 1 ";
			else if (tempp.getSymbol().equals(CustomEnum.inhibitionActivation)|| tempp.getSymbol().equals(CustomEnum.inhibitionInhibition) || tempp.getSymbol().equals(CustomEnum.activationInhibition)|| tempp.getSymbol().equals(CustomEnum.activationActivation))
				printx = printx + " 1 "+ " 1 ";
		}
		writeInFIle_NetworkFile(printx);

	}

	private static synchronized void writeInFIleMergedPaths() {

		Charset charset = StandardCharsets.UTF_8;
		try {
			List<Double> pavlueListString = new ArrayList<>();
			pavlueListString.clear();
			for (int sList = 0; sList < AllMergedpaths.size(); sList++) {
				Integer asd = sList;
				ArrayList<String> genesList = (ArrayList<String>) AllMergedpaths.get(sList);
				pavlueListString = new ArrayList<>();
				pavlueListString.clear();
				for (String string : genesList) {
					pavlueListString.add(genesKeyPValuesPair.get(string));
				}
				Double calculationResult = GenericFunctions.hartungFunction(pavlueListString);
				Files.write(Paths.get(Step3Output),
						(genesList.toString() + " " + calculationResult + "\n").getBytes(charset),
						StandardOpenOption.CREATE, StandardOpenOption.APPEND);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static List<String> getGenesList(Integer index) {
		try {
			return hashTableList.get(index);
		} catch (Exception e) {

			return new ArrayList<>();
		}

	}

	private static void setGenesList(List<String> list, Integer index) {
		hashTableList.put(index, list);
	}

	private static List<String> sortGeneInfoList(List<String> mainValueList) {
		if (mainValueList == null)
			return null;
		List<GenesInfo> tempList = new ArrayList<GenesInfo>();
		List<String> returnList = new ArrayList<>();
		for (String item : mainValueList) {
			GenesInfo temp = new GenesInfo();
			temp.setGeneName(item);
			try {
				temp.setpValue(genesKeyPValuesPair.get(item));
			} catch (Exception e) {
				temp.setpValue(0.0);
			}
			if (temp.getpValue() > 0 && temp.getpValue() < 1)
				tempList.add(temp);
		}

		Collections.sort(tempList);
		for (GenesInfo returnL : tempList) {
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
}