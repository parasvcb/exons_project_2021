import React, { Component } from "react";
const reactStringReplace = require("react-string-replace");

export class RenderExonsFull extends Component {
  constructor(props) {
    super(props);
    this.state = {
      showAaseq: true,
      reset: false,
      showSs: false,
      showDom: false,
      showAlert: false,
      showExob: false,
      showDis: false,
      showFasta: false
    };
    this.toggle_fasta = this.toggle_fasta.bind(this);
    this.toggle_aa = this.toggle_aa.bind(this);
    this.toggle_ss = this.toggle_ss.bind(this);
    this.toggle_dom = this.toggle_dom.bind(this);
    this.toggle_dis = this.toggle_dis.bind(this);

    var text = this.state.showAaseq
      ? "aaseq"
      : this.state.showDom
      ? "dom"
      : this.state.showDis
      ? "dis"
      : this.state.showSs
      ? "ss"
      : "fasta";

    this.props.handleUpdateInfRend(text);
  }

  toggle_ss() {
    this.setState({
      showAaseq: false,
      showDom: false,
      showSs: true,
      reset: true,
      showDis: false,
      showFasta: false
    });
    this.props.handleUpdateInfRend("ss");
  }
  toggle_dis() {
    this.setState({
      showAaseq: false,
      showDom: false,
      showSs: false,
      reset: true,
      showDis: true,
      showFasta: false
    });
    this.props.handleUpdateInfRend("dis");
  }

  toggle_dom() {
    this.setState({
      showAaseq: false,
      showDom: true,
      showSs: false,
      reset: true,
      showDis: false,
      showFasta: false
    });
    this.props.handleUpdateInfRend("dom");
  }

  toggle_aa() {
    this.setState({
      showAaseq: true,
      showDom: false,
      showSs: false,
      reset: false,
      showDis: false,
      showFasta: false
    });
    this.props.handleUpdateInfRend("aaseq");
  }
  toggle_fasta() {
    this.setState({
      showAaseq: false,
      showDom: false,
      showSs: false,
      reset: true,
      showDis: false,
      showFasta: true
    });
    console.log("clicked");
    this.props.handleUpdateInfRend("fasta");
  }

  render() {
    const exonHash = this.props.exOn;
    const trans = this.props.transwhole;
    const domset = this.props.domset;
    let aaVista=''
    trans.exonsIds.split(",").map(ex => (aaVista += exonHash[parseInt(ex)].aaseq));
    //console.log("aa",aaVista)

    // let domseqVista=''
    // trans.exonsIds.split(",").map(ex => (domseqVista += this.propConcate(trans.id,domset,exonHash[parseInt(ex)],"domseq")))

    let sseqVista=''
    trans.exonsIds.split(",").map(ex => (sseqVista += this.propConcate(trans.id,domset,exonHash[parseInt(ex)],"ssseq")))
    
    // let disseqVista=''
    // trans.exonsIds.split(",").map(ex => (disseqVista += this.propConcate(trans.id,domset,exonHash[parseInt(ex)],"disseq")))
    

    // console.log("domseq",domseqVista)
    console.log("ssseq",sseqVista)
    // console.log("disseq",disseqVista)
    
    // trans.exonsIds.split(",").map(ex => (domseqVista += exonHash[parseInt(ex)].exonsDom[0].domseq));

    // trans.exonsIds.split(",").map(ex => (domseqVista += exonHash[parseInt(ex)].exonsDom[0].domseq));
    
    //
    //console.log("TID",trans.tId)
    return (
      <>
        <span className="badge badge-info align-center">{trans.tId}</span>
        <div className="row">
          <div
            className="col"
            style={{
              flexDirection: "row",
              overflowX: "scroll"
            }}
          >
            <div
              className="d-flex flex-wrap"
              style={{ fontFamily: "monospace", wordWrap: "break-word" }}
            >
              {" "}
              {this.state.showFasta === true
                ? this.renderFasta(trans.tId, trans.exonsIds, exonHash)
                : null}
              {!this.state.showFasta
                ? trans.exonsIds.split(",").map(ex => (
                    <div key={parseInt(trans.id.toString() + ex)}>
                      <span
                        style={{ fontSize: 20, wordWrap: "break-word" }}
                        onMouseLeave={() => this.props.removeAlert()}
                        onMouseEnter={() =>
                          this.props.addAlert(
                            this.props.exonAlerts[exonHash[parseInt(ex)].exId]
                          )
                        }
                      >
                        {this.decideText(
                          trans.id,
                          domset,
                          exonHash[parseInt(ex)]
                        )}
                      </span>
                    </div>
                  ))
                : null}
            </div>
            <div className="container">
              {this.state.reset ? (
                <>
                  <button
                    onClick={this.toggle_aa}
                    className="m-2 btn btn-warning btn-md"
                  >
                    Reset to aa
                  </button>
                  <button className="m-2 btn btn-light btn-sm" disabled>
                    Fasta
                  </button>
                  <button className="m-2 btn btn-light btn-sm" disabled>
                    Show SS
                  </button>
                  <button className="m-2 btn btn-light btn-sm" disabled>
                    Show Dom
                  </button>
                  <button className="m-2 btn btn-light btn-sm" disabled>
                    Show Dis
                  </button>
                </>
              ) : (
                <>
                  <button className="m-2 btn btn-light btn-sm" disabled>
                    Reset to aa
                  </button>
                  <button
                    onClick={this.toggle_fasta}
                    className="m-2 btn btn-warning btn-md"
                  >
                    Fasta
                  </button>
                  <button
                    onClick={this.toggle_ss}
                    className="m-2 btn btn-warning btn-md"
                  >
                    Show SS
                  </button>
                  <button
                    onClick={this.toggle_dom}
                    className="m-2 btn btn-warning btn-md"
                  >
                    Show Dom
                  </button>
                  <button
                    onClick={this.toggle_dis}
                    className="m-2 btn btn-warning btn-md"
                  >
                    Show Dis
                  </button>
                </>
              )}
            </div>
          </div>
        </div>
      </>
    );
  }
  renderFasta(transid, transexonids, exonHash) {
    let fasta = "";
    let header = ">" + transid;
    let arr = [];
    arr.push(header);
    transexonids.split(",").map(ex => (fasta += exonHash[parseInt(ex)].aaseq));
    let uppfac = parseInt(fasta.length / 80);
    for (var i = 0; i <= uppfac; i++) {
      arr.push(fasta.slice(i * 80, (i + 1) * 80));
    }
    //console.log("fastaview,", arr);
    return (
      <span style={{ fontSize: 15, whiteSpace: "pre-line" }}>
        {arr.join("\n")}
      </span>
    );
  }
  transSisdomss(looper, trid, property) {
    //console.log("looper,", property);
    let tstate = true;
    if (looper.length < 2) {
      //console.log(looper.length === 0 ? "false" : looper[0][property]);
      return looper.length === 0 ? false : looper[0][property];
    } else {
      let patt = RegExp(`${trid}`, "g");
      for (var i = 0; i < looper.length; i++) {
        if (patt.exec(looper[i].list_trans_fk)) {
          //console.log(looper[i][property]);
          //console.log(property,looper[i][property])
          return looper[i][property];
        }
      }
    }

    if (tstate) {
      //console.log("false");
      return false;
    }
  }

  propConcate(transfk, domset, singleexon, visprop){
    //console.log("sing",singleexon)
    let transfkstr = transfk.toString();
    if (visprop=="domseq"){
    return this.transSisdomss(singleexon.exonsDom, transfk, visprop);
    }
    else if (visprop=="ssseq"){
      return this.transSisdomss(singleexon.exonsSS, transfk, visprop);
    }
    else {
      return this.transSisdomss(singleexon.exonsDis, transfk, visprop);
    }
  }

  decideText(transfk, domset, singleexon) {
    let transfkstr = transfk.toString();
    if (this.state.showAaseq) {
      return this.decorateText(
        singleexon.aaseq,
        singleexon.exId,
        [],
        "aaseq",
        []
      );
    }
    if (this.state.showSs) {
      let ss = this.transSisdomss(singleexon.exonsSS, transfk, "ssseq");
      return this.decorateText(ss, singleexon.exId, ["H", "E"], "ss", []);
    }
    if (this.state.showDis) {
      let dis = this.transSisdomss(singleexon.exonsDis, transfk, "disseq");
      return this.decorateText(dis, singleexon.exId, ["S"], "dis", []);
    }

    if (this.state.showDom) {
      let dom = this.transSisdomss(singleexon.exonsDom, transfk, "domseq");
      return this.decorateText(dom, singleexon.exId, [], "dom", domset);
    }
  }

  decorateText(textassum, id, vals, keyst, domset) {
    //let opaval = id[0] === "R" ? 0.25 : id.split(".")[2] === "A" ? 0.25 : id.split(".")[2] === "F" ? 0.5  : 0.7;
    //console.log("idexon", id, opaval);
    let replacedText = textassum.length > 0 ? textassum : false;
    //console.log(id, replacedText, vals, keyst);
    if (replacedText) {
      //console.log("yes replaed text");
      if (keyst === "ss") {
        vals.map(iter => {
          let color = iter === "H" ? "purple" : "yellow";
          let patt = new RegExp(`(${iter}+)`, "g");
          replacedText = reactStringReplace(replacedText, patt, (match, i) => (
            //console.log("ss", match, i, match + i, id),
            <span key={match + i} style={{ background: color, color: color }}>
              {match}
            </span>
          ));
        });
      } else {
        if (keyst === "dom") {
          domset.map(iter => {
            let patt = new RegExp(`(${iter.code}+)`, "g");
            replacedText = reactStringReplace(
              replacedText,
              patt,
              (match, i) => (
                <span
                  key={match + i}
                  style={{
                    background: iter.color,
                    color: iter.color
                  }}
                >
                  {match}
                </span>
              )
            );
          });
        } else {
          if (keyst === "dis") {
            let iter = "S";
            let patt = new RegExp(`(${iter}+)`, "g");
            replacedText = reactStringReplace(
              replacedText,
              patt,
              (match, i) => (
                //console.log("ss", match, i, match + i, id),
                <span
                  key={match + i}
                  style={{ background: "olive", color: "olive" }}
                >
                  {match}
                </span>
              )
            );
          }
        }
      }
    }
    if (replacedText) {
      //console.log("yesitis", replacedText, id);

      return (
        <span
          className="badge"
          style={{
            backgroundColor: "linen",
            color: "black",
            border: "1px solid black",
            marginRight: "10px"
          }}
        >
          {replacedText}
        </span>
      );
    } else {
      return (
        <span
          className="badge"
          style={{
            backgroundColor: "linen",
            color: "grey",
            width: "100px",
            border: "1px solid black",
            wordWrap: "break-word",
            marginRight: "10px"
          }}
        >
          UTR
        </span>
      );
    }
  }
}

export default RenderExonsFull;
