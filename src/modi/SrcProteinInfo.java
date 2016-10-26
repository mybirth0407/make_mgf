package modi;

import java.io.Serializable;

public class SrcProteinInfo implements Serializable {
  int srcProteinID;
  int startPos, endPos;

  public SrcProteinInfo(int srcProteinID, int startPos, int endPos) {
    this.srcProteinID = srcProteinID;
    this.startPos = startPos;
    this.endPos = endPos;
  }

  public int getSrcProteinID() {
    return srcProteinID;
  }

  public int getStartPos() {
    return startPos;
  }

  public int getEndPos() {
    return endPos;
  }

  public void setStartPos(int pos) {
    startPos = pos;
  }

  public String toString() {
    return "(SRC:\"" + srcProteinID + "\"" + startPos + "~" + endPos + ")";
  }
}
