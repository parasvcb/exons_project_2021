import React, { Component } from "react";

export class RenderExonsAln extends Component {
  // need to iterate the whole set together, not only single object

  listGiver(trans, exonhash) {
    let transexons = trans.exonsIds.split(",").map(function(x) {
      return parseInt(x, 10);
    });

    let isretention = false;
    const fillRange = (start, end) => {
      return Array(end - start + 1)
        .fill()
        .map((item, index) => start + index);
    };
    var arrayInt = fillRange(1, this.props.total);
    let arrayint_copy = arrayInt;

    Array.prototype.insert = function(index, item) {
      this.splice(index, 0, item);
    };
    var transexonmap = [];
    transexons.forEach(element => {
      if (exonhash[element].exId[0] !== "R") {
        transexonmap.push(parseInt(exonhash[element].exId.split(".")[3]));
      } else {
        isretention = true;
        transexonmap.push(exonhash[element].exId);
      }
    });
    var retention = [];
    if (isretention === true) {
      retention = transexonmap.filter(function(exon) {
        return exon[0] === "R";
      });
    }
    //transexonmap will hold the data of the current transcripts exons in a linear fashion
    //along iwth retention id where it is
    if (retention.length > 0) {
      var comMiniarr = [];
      retention.forEach(element => {
        var start = parseInt(element.split(":")[2].split(".")[3]);
        var end = parseInt(element.split(":")[4].split(".")[3]);
        var indextobeadded = arrayint_copy.indexOf(start);
        arrayint_copy.insert(indextobeadded, element);
        //adding retention exon to the position of all the exons
        //based on its start position in whole template positions
        //console.log("((((((((((index,arrold,newarr", indextobeadded, arrayint_copy);
        //console.log("startend", start, end);
        var miniarr = fillRange(start, end);
        comMiniarr = comMiniarr.concat(miniarr);
        //storing the position of spanned elements(exons) by retenion exon
        //console.log("miniarr0", miniarr);
      });
      //console.log("miniarr", comMiniarr);
      // if any retention exon, the above function will create a list of
      //spanning exon elements for reteion cases and concat them in comMiniarr
      //console.log("arrayintbefroe", arrayint_copy);
      arrayint_copy = arrayint_copy.filter(function(exon) {
        return comMiniarr.indexOf(exon) === -1;
      });
      //elemnts thar are spanned by retenion cases ae removed from arrayint_copy
    }
    //console.log("arrayintAfter", arrayint_copy);

    // the above block will filter arrayint_copy which is all exons for the gene
    // in this if any retenion spanned exons are there, from comMiniarr
    // then they are automatically removed, as we will filter out all the
    // elements of arryINt which are not preset in the coMiniarr.
    var finalwriteup = [];
    for (var i = 0; i < arrayint_copy.length; i++) {
      var tempob = {};
      for (var j = 0; j < transexons.length; j++) {
        if (typeof arrayint_copy[i] === "number") {
          if (
            exonhash[transexons[j]].exId[0] !== "R" &&
            parseInt(exonhash[transexons[j]].exId.split(".")[3]) ===
              arrayint_copy[i]
          ) {
            tempob = exonhash[transexons[j]];
            break;
          }
        } else {
          if (exonhash[transexons[j]].exId === arrayint_copy[i]) {
            tempob = exonhash[transexons[j]];
            break;
          }
        }
      }
      if (Object.keys(tempob).length === 0 && tempob.constructor === Object) {
        tempob = { exId: arrayint_copy[i] };
      }
      finalwriteup.push(tempob);
    }
    return finalwriteup;
  }

  // i will have a list with elelmnt t be renederd with empty number
  // render with name ans rest for empty candidadtes,
  // render the dotted output
  // for ncb cases reder the boundaries dotted
  //

  render() {
    //console.log("props->", this.props);

    const exonhash = this.props.exOn;
    var exonset = [];
    this.props.arraylis.map(trans => {
      let arrayout = this.listGiver(trans, exonhash);
      let temparr = { array: arrayout, tid: trans.tId };
      exonset.push(temparr);
      //console.log("trans output", trans, "out->", exonset);
    });

    //retention array will keep a count on the retention cases if any
    //in the transexonmap

    //console.log("transexonsmap", transexonmap, arrayInt, finalwriteup);
    let ext = "";
    return (
      <div
        className="container-fluid"
        style={{ overflowX: "auto", padding: "1rem" }}
      >
        {exonset.map(tran => (
          <div key={tran.tid} style={{ display: "flex", padding: "0.5rem" }}>
            <span>
              <span className="badge badge-success" style={{ width: "150px" }}>
                {tran.tid}->
              </span>
            </span>
            {tran.array.map(ex => (
              <span key={ex.exId}>
                <span
                  className="badge "
                  style={this.alnStyle(ex.exId, this.props.exonCodes)}
                >
                  {ex.exId}
                  {this.ncbStyler(ex.exId)}
                </span>
              </span>
            ))}
          </div>
        ))}
      </div>
    );
  }
  ncbStyler(exid) {
    let typeid = typeof exid;
    let ob = null;
    if (typeid === "string") {
      if (exid[0] === "R") {
        let elem = exid.split(":");
        let first = elem[2].split(".")[4];
        let second = elem[4].split(".")[4];
        let elemneeded = first + second;
        switch (elemneeded) {
          case "n0":
            ob = (
              <span className="badge badge-light float-left border">{"<"}</span>
            );
            break;
          case "0c":
            ob = (
              <span className="badge badge-light float-right border">
                {">"}
              </span>
            );
            break;

          case "nc":
            ob = (
              <>
                <span className="badge badge-light float-left border">
                  {"<"}
                </span>
                <span className="badge badge-light float-right border">
                  {">"}
                </span>
              </>
            );
            break;
          default:
            ob = null;
        }
      } else {
        let ext = exid.split(".")[4];
        switch (ext) {
          case "0":
            ob = null;
            break;
          case "n":
            ob = (
              <span className="badge badge-light float-left border">{"<"}</span>
            );
            break;
          case "c":
            ob = (
              <span className="badge badge-light float-right border">
                {">"}
              </span>
            );
            break;

          case "b":
            ob = (
              <>
                <span className="badge badge-light float-left border">
                  {"<"}
                </span>
                <span className="badge badge-light float-right border">
                  {">"}
                </span>
              </>
            );
            break;
          default:
            ob = null;
        }
      }
    }
    return ob;
  }
  alnStyle(exid, codes) {
    let wid = "100px";
    if (typeof exid === "number") {
      let color = "";
      //console.log("color", exid, codes, codes[exid]);
      switch (codes[exid]) {
        case "T":
          color = "darkgreen";
          break;
        case "D":
          color = "deeppink";
          break;
        case "U":
          color = "goldenrod";
          break;
        case "M":
          color = "lightsalmon";
          break;
        default:
          color = "white";
      }
      return {
        width: wid,
        marginLeft: "50px",
        marginRight: "50px",
        backgroundColor: color,
        color: "white"
      };
    } else {
      let color = "";
      let testval = "";
      if (exid[0] === "R") {
        let ex = exid.split(":");
        let st = parseInt(ex[2].split(".")[3]);
        let ed = parseInt(ex[4].split(".")[3]);
        wid = ((ed - st + 2) * 100).toString() + "px";
        testval = exid.split(":")[1];
      }
      if (testval.length === 0) {
        testval = exid.split(".")[1];
      }

      switch (testval) {
        case "-2":
          color = "goldenrod";
          break;

        case "1":
          color = "darkgreen";
          break;
        case "-1":
          color = "goldenrod";
          break;
        case "0":
          color = "cornsilk";
          break;
        case "2":
          color = "darkolivegreen";
          break;
        case "5":
          color = "darkseagreen";
          break;

        case "4":
          color = "seagreen";
          break;
        case "3":
          color = "forestgreen";
          break;
        case "6":
          color = "mediumseagreen";
          break;
        case "7":
          color = "limegreen";
          break;
        case "8":
          color = "palegreen";
          break;
        case "9":
          color = "lightgreen";
          break;

        default:
          color = "white";
      }
      return {
        width: wid,
        marginLeft: "50px",
        marginRight: "50px",
        backgroundColor: color,
        color: "white"
      };
    }
  }
}

export default RenderExonsAln;
