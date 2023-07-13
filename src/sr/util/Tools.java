package sr.util;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import sr.evolution.tree.SRTree;

import java.util.*;

public class Tools {

	public final static double globalPrecisionThreshold = 1e-10;

	public static String getFirstLeafID(Node node) {
		if (node.isLeaf()) {
			return node.getID();
		} else {
			List<Node> leafs = node.getAllLeafNodes();
			Collections.sort(leafs, Comparator.comparingDouble(Node::getHeight).reversed());
			return leafs.get(0).getID();
		}
	}

	public static String removeLastSubstring(String split, String s) {
		String[] substrings = s.split(split);
		String[] substringsExceptLast = Arrays.copyOfRange(substrings, 0, substrings.length - 1);
		return String.join(split, substringsExceptLast);
	}

	/**
	 * Orientate node children depending on stored metadata.
	 *
	 * @param subtreeRootNr the none number
	 */
	public static void orientateNodeChildren(int subtreeRootNr, Tree tree) {
		Node subTreeRoot = tree.getNode(subtreeRootNr);
		if (!subTreeRoot.isLeaf()) {
			if(subTreeRoot.getChild(0).metaDataString == null)
				return;
			if ((!subTreeRoot.isFake() && !subTreeRoot.getChild(0).metaDataString.contains("orientation=ancestor"))
					|| (subTreeRoot.isFake() && subTreeRoot.getChild(1).getHeight() != subTreeRoot.getHeight())) {
				Node left = subTreeRoot.getChild(1);
				Node right = subTreeRoot.getChild(0);

				subTreeRoot.removeAllChildren(false);

//                subTreeRoot.setLeft(left);
//                subTreeRoot.setRight(right);
				subTreeRoot.addChild(left);
				subTreeRoot.addChild(right);
			}

			orientateNodeChildren(subTreeRoot.getChild(0).getNr(), tree);
			orientateNodeChildren(subTreeRoot.getChild(1).getNr(), tree);
		}
	}

	public static boolean equalHeightWithPrecision(Node n1, Node n2) {
		if (Math.abs(n1.getHeight() - n2.getHeight()) <= globalPrecisionThreshold)
			return true;
		return false;
	}

	public static Set<String> taxaNamesToSpeciesNames(Set<String> taxonNames, String splitCharacter){
		Set<String> speciesNames = new HashSet<>();
		for(String taxonName : taxonNames){
			speciesNames.add(removeLastSubstring(splitCharacter, taxonName));
		}
		List<String> list = new ArrayList<>(speciesNames);
		Collections.sort(list); // Sort the list alphabetically

		Set<String> sortedSpeciesNames = new LinkedHashSet<>(list); // Convert the sorted list back to a set

		return sortedSpeciesNames;
	}

	public static Integer getSpeciesNumber(Set<String> speciesNames, String id ){
		Iterator<String> it = speciesNames.iterator();
		int i = 0;
		while (it.hasNext()){
			if(it.next().contains(id))
				return i;
			i++;
		}
		return -1;
	}

	public static double[] getSeparatingLengthAndNodeCount(Node parent, Node child){
		double length = 0.;
		double nodeCount =0.;
		length += child.getLength();
		Node tmpParent = child.getParent();
		while (parent.getNr() != tmpParent.getNr()){
			nodeCount++;
			child = tmpParent;
			length += child.getLength();
			tmpParent = tmpParent.getParent();
		}
		return new double[]{length, nodeCount};
	}

	private static void increaseLastDouble(List<Double> list, double value) {
		list.set(list.size() - 1, list.get(list.size() - 1) + value);
	}

}
