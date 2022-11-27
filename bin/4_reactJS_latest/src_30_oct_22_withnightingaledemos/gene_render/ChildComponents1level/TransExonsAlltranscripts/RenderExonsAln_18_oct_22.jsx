import React, { Component } from "react";

export class RenderExonsAln extends Component {
  // need to iterate the whole set together, not only single object

  listGiver(trans, exonhash) {
        let transexons = trans.exonsIds.split(",").map(function(x) {
          return parseInt(x, 10);
        });
        // split the oferign keys and get there strings in the list

        let isretention = false;

        const fillRange = (start, end) => {
          return Array(end - start + 1)
            .fill()
            .map((item, index) => start + index);
        };
        // justa function definition

        var arrayInt = fillRange(1, this.props.total);
        //total is final exon in geneset

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
      //console.log("trans output",
      // split the oferign keys and get there strings in the list trans, "out->", exonset);
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
          <div key={tran.tid} style={{ display: "flex", padding: "0.5rem"}}>
            <span>
              <span className="badge bg-success" style={{ width: "150px" }}>
                <span>
                {tran.tid}
                </span>
              </span>
            </span>
            {tran.array.map(ex => (
              <span key={ex.exId}>
                  <span 
                    className="badge"
                    style={this.alnStyle(ex.exId, this.props.exonCodes)}
                  >
                        {(this.ncbStylerLeft(ex.exId) !== null & this.ncbStylerRight(ex.exId)!=null)
                         ? //Both N and C
                              <div style={{
                                  fontSize:"12px",           
                                  fontWeight: "bold",
                                  borderLeft:"5px solid #D06224",
                                  padding: "0px 2px 0px 2px", 
                                  borderRight:"5px solid #D06224"}}>
                                  {ex.exId}
                              </div>
                              // only N case below
                         : this.ncbStylerLeft(ex.exId) !== null ?
                              <span style={{
                                fontSize:"12px",
                                fontWeight: "bold",
                                padding: "0px 5px 0px 5px", 
                                borderLeft:"5px solid #D06224", 
                                borderRight:"5px solid transparent"}}>
                                {ex.exId}
                              </span>
                              // only C case below
                         : this.ncbStylerRight(ex.exId) !== null ?
                              <span style={{
                                fontSize:"12px",
                                fontWeight: "bold",
                                padding: "0px 5px 0px 5px", 
                                borderLeft:"5px solid transparent", 
                                borderRight:"5px solid #D06224"}}>
                                {ex.exId}
                              </span>
                        : <span style={{
                          fontSize:"15px",
                          fontWeight: "lighter"}}>
                        {ex.exId}</span>}   
              </span>
              </span>
              
            ))}
          </div>
        ))}
      </div>
    );
  }
  ncbStylerLeft(exid) {
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
                <i className= "bi bi-skip-start-fill m-0 p-0"   style = {{ fontSize:"15px", color: "#D06224"}}>   </i>
              
            );
            break;

          case "nc":
            ob = (
                  <i className= "bi bi-skip-start-fill m-0 p-0"   style = {{ fontSize:"15px", color: "#D06224"}}>   </i>
              
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
                <i className= "bi bi-skip-start-fill m-0 p-0"   style = {{ fontSize:"15px", color: "#D06224"}}>   </i>
            );
            break;
          case "b":
            ob = (
                <i className= "bi bi-skip-start-fill m-0 p-0"   style = {{ fontSize:"15px", color: "#D06224"}}>   </i>
            );
            break;
          default:
            ob = null;
        }
      }
    }
    return ob;
  }
  ncbStylerRight(exid) {
    let typeid = typeof exid;
    let ob = null;
    if (typeid === "string") {
      if (exid[0] === "R") {
        let elem = exid.split(":");
        let first = elem[2].split(".")[4];
        let second = elem[4].split(".")[4];
        let elemneeded = first + second;
        switch (elemneeded) {
          case "0c":
            ob = (
                <i className= "bi bi-skip-end-fill m-0 p-0"   style = {{ fontSize:"15px", color: "white"}}>   </i>
            );
            break;

          case "nc":
            ob = (
                  <i className= "bi bi-skip-end-fill m-0 p-0"   style = {{fontSize:"15px", color: "white"}}>   </i>
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
          case "c":
            ob = (
                <i className= "bi bi-skip-end-fill m-0 p-0"   style = {{fontSize:"15px", color: "white"}}>   </i>
            );
            break;

          case "b":
            ob = (
                  <i className= "bi bi-skip-end-fill"   style = {{fontSize:"15px", color: "white"}}>   </i>
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
        backgroundColor: "white",
        border: "dotted 3px " + "black",
        color: "black",
      };
    } else {
      let color = "white";
      let background='';
      let testval = "";
      if (exid[0] === "R") {
        let ex = exid.split(":");
        let st = parseInt(ex[2].split(".")[3]);
        let ed = parseInt(ex[4].split(".")[3]);
        // wid = ((ed - st + 2) * 100).toString() + "px";
        wid = (((ed - st + 1)+(ed-st)) * 100).toString() + "px";
        testval = exid.split(":")[1];
      }
      if (testval.length === 0) {
        testval = exid.split(".")[1];
      }


      switch (testval) {
        case "-2":
          background = "#cdcdcd";
          color ="black"
          break;
        case "-1":
          background = "goldenrod";
          break;
        case "0":
          background = "cornsilk";
          color = "black";
          break;
        case "1":
          // background = "darkgreen";
          background = "#CFFF8D"
          color = "black";
          break;
        case "2":
          background = "#A8E890";
          color = "black";
          break;
        case "3":
          background = "#9CFF2E";
          color = "black"
          break;
        case "4":
          background = "#38E54D";
          color = "black";
          break;
        case "5":
          background = "#749F82";
          color = "black";
          break;  
        case "6":
          background = "#3D8361";
          color = "white"
          break;
        case "7":
          background = "#1C6758";
          color = "white"
          break;
        case "8":
          background = "#425F57";
          color = "white"
          break;
        case "9":
          background = "#483838";
          color = "white"
          break;
        default:
          background = "#000000";
          color = "white"
      }
      return {
        width: wid,
        marginLeft: "50px",
        marginRight: "50px",
        backgroundColor: background,
        color: color
      };
    }
  }
}
/*
// {this.ncbStylerLeft(ex.exId) !== null ? <span  style = {{borderColor:"#D06224", borderLeft:"5px solid #D06224"}}>{this.ncbStylerLeft(ex.exId)}</span> : <span  style = {{borderColor:"#D06224", borderLeft:"5px solid transparent"}}>{this.ncbStylerLeft(ex.exId)}</span>}
// <span >{ex.exId}</span>
// {this.ncbStylerRight(ex.exId) !== null ? <span className="border-end border-5 float-right">{this.ncbStylerRight(ex.exId)}</span> : <span  style = {{borderRight:"5px solid transparent"}}>{this.ncbStylerLeft(ex.exId)}</span>}
// </span>

{this.ncbStylerLeft(ex.exId) !== null ? <span  style = {{borderColor:"#D06224", borderLeft:"5px solid #D06224"}}>{this.ncbStylerLeft(ex.exId)}</span> : <span  style = {{borderColor:"#D06224", borderLeft:"5px solid transparent"}}>{this.ncbStylerLeft(ex.exId)}</span>}
                   <span style={{fontSize:"10px"}} >{ex.exId}</span>
                   {this.ncbStylerRight(ex.exId) !== null ? <span className="border-end border-5 float-right">{this.ncbStylerRight(ex.exId)}</span> : <span  style = {{borderRight:"5px solid transparent"}}>{this.ncbStylerLeft(ex.exId)}</span>}
*/
export default RenderExonsAln;


