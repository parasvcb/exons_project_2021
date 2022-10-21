import React, { Component } from "react";
let jsdiff = require("diff");

export class SeqExonMatching extends Component {
  constructor(props) {
    super(props);
    this.state = { alert: false };
    this.transcriptOrderedArray = this.transcriptOrderedArray.bind(this);
    this.firstpassArrfunc = this.firstpassArrfunc.bind(this);
    this.finalExonsList = this.finalExonsList.bind(this);
  }

  addAlert(exms) {
    this.setState({
      alert: exms
    });
  }
  removeAlert() {
    this.setState({
      alert: false
    });
  }

  transcriptOrderedArray(transprops) {
    //transprops here is the exongenes
    let allaa = {};
    console.log("transpropsin()", transprops);
    transprops.forEach(exons => {
      if (exons.length > 0) {
        let aaseq = "";
        let position = "";
        let exid = exons.exId;
        if (exons.exId[0] === "R") {
          position = exons.exId.split(":")[2].split(".")[3];
        } else {
          position = exons.exId.split(".")[3];
        }
        aaseq = exons.aaseq;

        if (allaa.hasOwnProperty(position)) {
          //console.log("yes in if", position, allaa);
          if (!allaa[position].hasOwnProperty(exid)) {
            allaa[position][exid] = aaseq;
            //console.log("yes in if2", position, allaa);
          }
        } else {
          //console.log("yes in else", position, allaa);
          allaa[position] = {};
          allaa[position][exid] = aaseq;
          //console.log("yes in if2", position, allaa);
        }
      }
    });
    let arr = [];
    for (let key in allaa) {
      let testob = {};
      testob.exonslist = [];
      testob.exindex = parseInt(key);
      for (let exon in allaa[key]) {
        let t_ob = {};
        t_ob.exonid = exon;
        t_ob.aaseq = allaa[key][exon];
        testob.exonslist.push(t_ob);
      }
      arr.push(testob);
    }
    //console.log("allaa", allaa);
    console.log("allaa_modified", arr);
    arr.sort(function(a, b) {
      return a.exindex - b.exindex;
      //sorted in ascending order
    });
    /*
    arr is the array having objects sorted in the asceding order of exon ids
    each object
    {
      exonslist:[],
      exindex:int(0)
    } 
    has exonslist which further has object_level2 with atributes 
    {
      exonid:str,
      aaseq:str,
    }
    */
    return arr;
  }

  firstpassArrfunc(transOrdered, userfasta) {
    let sequence = userfasta;
    let matrix = [];
    let anchor = transOrdered[0].exindex;
    transOrdered.forEach(position => {
      let dummysort = [];
      position.exonslist.forEach(exon => {
        let patt = new RegExp(`${exon.aaseq}`);
        let match = patt.exec(sequence);
        if (match) {
          //console.log("matched", match);
          let testob = {
            exonid: exon.exonid,
            len: exon.aaseq.length,
            index: match.index,
            aaseq: exon.aaseq
          };
          dummysort.push(testob);
        }
      });
      //need a else equivalent of the following, (29/7/19)
      //not needed, as empty will be the one exon group,  with no match in userfasta
      if (dummysort.length > 0) {
        dummysort.sort(function(a, b) {
          return b.len - a.len;
          //sorted in descending order of matched length
        });
        //console.log("InDummysort");
        let elemneed = dummysort[0];
        let insertOb = {};
        let unmataa = false;
        let obanchor = false;
        insertOb.exonid = elemneed.exonid;
        insertOb.exindex = position.exindex;
        insertOb.aaseq = elemneed.aaseq;

        if (elemneed.index !== 0) {
          unmataa = sequence.slice(0, elemneed.index);
          obanchor = [anchor, position.exindex];
        }
        anchor = position.exindex;
        insertOb.anchor = obanchor;
        insertOb.unmataa = unmataa;

        matrix.push(insertOb);
        //console.log("sequence before", sequence);

        sequence = sequence.slice(
          elemneed.index + elemneed.aaseq.length,
          sequence.length
        );
        //console.log("sequence after", sequence);
      }
    });
    if (sequence.length !== 0) {
      let testob = {};
      testob.anchor = [anchor, transOrdered[transOrdered.length - 1].exindex];
      testob.exindex = false;
      testob.aaseq = false;
      testob.unmataa = sequence;
      testob.exonid = false;
      matrix.push(testob);
      /*
      each element in the matrix will have a list whwe each list element will correspond to that of the object
      with follwing attributes
      {
        anchor:[from,to], range will exclude continuous stretche of matches (complete/partial) and will only include
        unmatched exons just after exons with matched(complete/partial) sequence
        exindex: numeric position of exon, (null case in last unmatched stretch, no match)
        aaseq:aminoacd sequnec, false(null case in last unmatched stretch, no match)
        unmataa: unmatched aminoacid sequence , flase(in perfect matches)
        exonid: exonid, (false,(null case in last unmatched stretch, no match))
      }
      */
    }

    return matrix;
  }
  finalExonsList(matrix, transOrdered) {
    //preliminary check if only one matrix elemnet is there
    const fillRange = (start, end) => {
      return Array(end - start + 1)
        .fill()
        .map((item, index) => start + index);
    };
    let final_list = [];
    matrix.forEach(exon => {
      let testob = {};
      if (!exon.unmataa) {
        testob.aaseq = exon.aaseq;
        testob.unmataa = false;
        testob.value = false;
        testob.exonid = exon.exonid;
        testob.status = 1;
        testob.message = "original";
        testob.len = exon.aaseq.length;
        final_list.push(testob);
      } else {
        let dummy_list = [];
        let arrayInt = fillRange(exon.anchor[0], exon.anchor[1]);
        if (arrayInt.length > 0) {
          let testob_default = {};
          testob.aaseq = false;
          testob.exonid = false;
          testob.unmataa = exon.unmataa;
          testob.value = false;
          testob.status = -1;
          testob.message =
            "Not matching to exon most probably it is, either C extension of former exon or N extension of latter with less chances of new insertion";

          //add a conditon here to watch, first and last cases
          testob.len = exon.unmataa.length;
          final_list.push(testob);
        } else {
          arrayInt.forEach(exonum => {
            if (transOrdered.hasOwnProperty(exonum)) {
              transOrdered[exonum].exonslist.forEach(exon2 => {
                let aaexon = exon2.aaseq;
                let matchchar = jsdiff.diffChars(aaexon, exon.unmataa);
                matchchar.forEach(chararr => {
                  if (
                    !(
                      chararr.hasOwnProperty("added") ||
                      chararr.hasOwnProperty("removed")
                    ) &&
                    chararr.count > 5
                  ) {
                    let pattern = new RegExp(`${chararr.value}`);
                    let match = pattern.exec(exon.unmataa);
                    let status = false;
                    let message = "";
                    if (match) {
                      message =
                        match.index === 0
                          ? "probable C-ter change"
                          : match.index === exon.unmataa.length - chararr.count
                          ? "probable N-ter change"
                          : "probable B case or new insertion";
                    }
                    testob.aaseq = exon2.aaseq;
                    testob.unmataa = exon.unmataa;
                    //here unmataa will be from left userdefonedsequnec's unmataa
                    testob.value = match.value;
                    //userdefined sequence value and unmat aa common set
                    testob.exonid = exon.exonid;
                    testob.status = 0;
                    testob.message = message;
                    testob.len = chararr.count;

                    dummy_list.push(testob);
                  }
                });
              });
            }
          });
          if (dummy_list.length > 0) {
            dummy_list.sort(function(a, b) {
              return b.len - a.len;
            });
            final_list.push(dummy_list);
          } else {
            let testob_default = {};
            testob.aaseq = false;
            testob.exonid = false;
            testob.unmataa = exon.unmataa;
            testob.value = false;
            testob.status = -1;
            testob.message =
              "no match found for this very region, likely new insertion or extension";
            testob.len = exon.unmataa.length;
            final_list.push(testob);
            //unmataa and nothing left to be matched in eoxns, with partial complete match, likely n extension
          }
        }
        //following condition to push only genuine cases after unmat aa and exclude empty push in last null case
        if (exon.aaseq) {
          let testob_ori = {};
          testob_ori.aaseq = exon.aaseq;
          testob_ori.exonid = exon.exonid;
          testob_ori.status = 1;
          testob_ori.unmataa = false;
          testob.value = false;
          testob_ori.message = "original";
          testob_ori.len = exon.aaseq.length;
          final_list.push(testob_ori);
          //after unmataa, push the original matched aaseq sequnec
        }
      }
    });
    return final_list;
    /*
    final list have in the end a list of exons and representative sequences
    with each object having follwing atributes
    {
      aaseq:sequnece of exon if status:true
      render if stauts is 1
      unmataa: render if status:0, this will be the unmatched component, with possible matches
      render it again if status:-1, NO mTACH forund
      value: show the match of exon aa to that of the left right exon
      exonid: given if stuats is 1 and 0, but no for -1
      status:true(perfectmatch),false()


    }
    */
  }

  transSeq(transobs, fasta, exonhash) {
    fasta = fasta.trim();
    let hasStru = {};
    //console.log("exonhas", exonhash);
    transobs.forEach(position => {
      let dummyseq = "";
      position.exonsIds.split(",").forEach(exon => {
        dummyseq += exonhash[parseInt(exon)].aaseq;
      });
      dummyseq = dummyseq.trim();
      //console.log("position tid", position);
      if (hasStru.hasOwnProperty(dummyseq)) {
        hasStru[dummyseq].push(position.tId);
      } else {
        hasStru[dummyseq] = [position.tId];
      }
    });
    //console.log("hasstructure", hasStru);
    console.log("fasta", fasta);
    console.log("check", hasStru[fasta]);
    if (hasStru.hasOwnProperty(fasta)) {
      return hasStru[fasta];
    } else {
      return false;
    }
  }

  render() {
    console.log("thispropstranscripts", this.props);
    //scan and fix the pre transcript ssequnece matches
    let transcripts_seq = this.transSeq(
      this.props.transcripts,
      this.props.fasta,
      this.props.exoncodes
    );
    console.log("transcriptsseq", transcripts_seq, typeof transcripts_seq);
    if (typeof transcripts_seq === "object") {
      let val = transcripts_seq.length > 1 ? "transcripts" : "transript";
      return (
        <>
          <h4>
            Your entered sequence, has identical protein product to following{" "}
            {val}:
          </h4>
          {transcripts_seq.map(ex => (
            <span
              key={ex}
              className="badge badge-light m-2"
              style={{ fontSize: 20 }}
            >
              {ex}
            </span>
          ))}
        </>
      );
    } else {
      let requiredtrans = this.transcriptOrderedArray(this.props.exOn);
      let matrix = this.firstpassArrfunc(requiredtrans, this.props.fasta);
      console.log("maths", this.props);
      console.log("required", requiredtrans);
      console.log("matrix1", matrix);
      // let match = jsdiff.diffChars("paras verma paras", "paras unknonw paras");
      // console.log(match);
      let exonsinlista = this.finalExonsList(matrix, requiredtrans);
      console.log("final_list", exonsinlista);
      return (
        <div className="row">
          <div
            className="col-sm-9 border"
            style={{
              flexDirection: "row",
              overflowX: "scroll"
            }}
          >
            <div
              className="d-flex flex-wrap"
              style={{ fontFamily: "monospace" }}
            >
              {exonsinlista.map(
                ex => (
                  console.log("tostring", ex.exonid, ex.unmataa, ex.aaseq),
                  (
                    <div
                      key={
                        ex.exonid.toString() +
                        ex.unmataa.toString() +
                        ex.aaseq.toString()
                      }
                    >
                      <span
                        onMouseLeave={() => this.removeAlert()}
                        onMouseEnter={() =>
                          this.addAlert(
                            ex.status === 1
                              ? this.props.exonAlerts[ex.exonid]
                              : ex.message
                          )
                        }
                      >
                        {this.decideText(ex)}
                      </span>
                    </div>
                  )
                )
              )}
            </div>
          </div>
          <div className="col-sm-3 border">
            <div className="d-flex flex-column sticky-top">
              <div className="p-2 bg-light">{this.showInfFunc()}</div>
            </div>
          </div>
        </div>
      );
    }
  }
  showInfFunc(val) {
    return this.state.alert ? (
      typeof this.state.alert === "string" ? (
        <p>
          <span
            className="badge badge-pill mr-3"
            style={{ whiteSpace: "pre-line" }}
          >
            {this.state.alert}
          </span>
        </p>
      ) : (
        this.state.alert
      )
    ) : null;
  }

  decideText(ex) {
    let sta = ex.status;
    let text = sta === 1 ? ex.aaseq : ex.unmataa;
    let bgc = sta === 1 ? "grey" : sta === 0 ? "silver" : "ivory";
    let col = sta === 1 || sta === 0 ? "white" : "black";
    return (
      <span className="badge m-2" style={{ backgroundColor: bgc, color: col }}>
        {text}
      </span>
    );
  }
}

export default SeqExonMatching;
