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

}
