import React, { Component } from "react";
import { decideText, propConcate, renderFasta } from "./ModuleREF_vista";

// const reactStringReplace = require("react-string-replace");
// This is updated file and should be worked on
// 
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
    // trans.exonsIds.split(",").map(ex => (domseqVista += propConcate(trans.id,domset,exonHash[parseInt(ex)],"domseq")))

    let sseqVista=''
    trans.exonsIds.split(",").map(ex => (sseqVista += propConcate(trans.id,domset,exonHash[parseInt(ex)],"ssseq")))
    
    // let disseqVista=''
    // trans.exonsIds.split(",").map(ex => (disseqVista += propConcate(trans.id,domset,exonHash[parseInt(ex)],"disseq")))
    // console.log("domseq",domseqVista)

    //console.log("ssseq",sseqVista)
    // console.log("disseq",disseqVista)
    // trans.exonsIds.split(",").map(ex => (domseqVista += exonHash[parseInt(ex)].exonsDom[0].domseq));
    // trans.exonsIds.split(",").map(ex => (domseqVista += exonHash[parseInt(ex)].exonsDom[0].domseq));
    
    //
    //console.log("TID",trans.tId)
    return (
      <div style= {{border: "double 1px lightgray", padding:"10px 5px 5px 10px"}}>
        <div className="h5">
        <span>Isoform:</span>
        <span className="badge align-center mx-2 my-3 h5"
              style = {{color: "#00303f", background: "#e9e9e9", borderLeft: "5px solid #00303f"}}
        >{trans.tId}</span>
        </div>
        <div className="row">
          <div className="col"
               style={{
                  flexDirection: "row",
                  overflowX: "scroll"
                }}
          >
            <div
              className="d-flex flex-wrap"
              style={{ fontFamily: "monospace", wordWrap: "break-word"}}
            >
              {" "}
              {this.state.showFasta === true
                ? renderFasta(trans.tId, trans.exonsIds, exonHash)
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
                        {decideText(
                          this.state,
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
                    className="m-2 btn btn-md"
                    style = {{background: "#b37d4e", color: "white"}}
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
                    className="mx-2 my-0 px-2 py-1 btn btn-md"
                    style= {{background:"#a7d2cb"}}
                  >
                    Fasta
                  </button>
                  <button
                    onClick={this.toggle_ss}
                    className="m-2 btn px-2 py-1 btn-md"
                    style= {{background:"#a7d2cb"}}
                  >
                    Show SS
                  </button>
                  <button
                    onClick={this.toggle_dom}
                    className="m-2 px-2 py-1 btn btn-md"
                    style= {{background:"#a7d2cb"}} //caebf2
                  >
                    Show Dom
                  </button>
                  <button
                    onClick={this.toggle_dis}
                    className="m-2 px-2 py-1 btn btn-md"
                    style= {{background:"#a7d2cb"}}
                  >
                    Show Dis
                  </button>
                </>
              )}
            </div>
          </div>
        </div>
      </div>
    );
  }

}

export default RenderExonsFull;
