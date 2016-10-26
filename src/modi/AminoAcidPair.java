package modi;

import java.util.ArrayList;
import java.util.Collections;

public class AminoAcidPair implements Comparable<AminoAcidPair> {
  static {
    AminoAcidPairList = new ArrayList<AminoAcidPair>();
    makeAAPairTable();
  }

  private AminoAcid firstAA;
  private AminoAcid secondAA;
  private double mass;

  public AminoAcidPair(AminoAcid firstAA, AminoAcid secondAA, double mass) {
    this.firstAA = firstAA;
    this.secondAA = secondAA;
    this.mass = mass;
  }

  public AminoAcid getFirstAA() {
    return firstAA;
  }

  public AminoAcid getSecondAA() {
    return secondAA;
  }

  public static ArrayList<AminoAcidPair> getCode(double mass, double tolerance) {
    ArrayList<AminoAcidPair> result = new ArrayList<AminoAcidPair>();

    for (int i = 0; i < AminoAcidPairList.size(); i++) {
      if (AminoAcidPairList == null)
        continue;
      else if (AminoAcidPairList.get(0).mass < mass - tolerance)
        break;
      else if (AminoAcidPairList.get(19).mass > mass + tolerance)
        break;
      else
        result.add(AminoAcidPairList.get(i));
    }

    return result;
  }

  public int compareTo(AminoAcidPair tag) {
    if (this.mass > tag.mass)
      return 1;
    else if (this.mass == tag.mass)
      return 0;
    else
      return -1;
  }

  public String toString() {
    return "" + firstAA + secondAA + " : " + mass;
  }

  private static final ArrayList<AminoAcidPair> AminoAcidPairList;

  private static void makeAAPairTable() {
    AminoAcid firstAA, secondAA;

    firstAA = AminoAcid.getAminoAcid((char) ('P'));
    for (int j = 0; j < 'Z' - 'A'; j++) {
      secondAA = AminoAcid.getAminoAcid((char) ('A' + j));
      if (secondAA == null)
        continue;
      AminoAcidPairList.add(new AminoAcidPair(firstAA, secondAA, firstAA.getMonoMass() + secondAA.getMonoMass()));
    }

    Collections.sort(AminoAcidPairList);
  }
}
