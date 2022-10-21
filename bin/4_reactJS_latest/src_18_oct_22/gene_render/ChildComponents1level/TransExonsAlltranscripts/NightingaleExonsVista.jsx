import React, { Component, Fragment, useState, useEffect, useCallback } from "react";
import ProtvistaManager from "protvista-manager";
import ProtvistaTrack from "protvista-track";
import ProtvistaNavigation from "protvista-navigation";
import ProtvistaSequence from "protvista-sequence";
import ProtvistaInterproTrack from "protvista-interpro-track";
import ProtvistaColouredSequence from "protvista-coloured-sequence";
import loadWebComponent from "./utils/load-web-component";
//import sequence from "../mocks/sequence.json";
// import { dataIPR, signatures } from "../mocks/interpro";
// import secondaryStructureData from "../mocks/interpro-secondary-structure.json";
import ProtvistaSaver from "protvista-saver";
import ProtvistaOverlay from "protvista-overlay";
import ProtvistaZoomTool from "protvista-zoom-tool";
import ProtvistaDatatable from "protvista-datatable";
//import { dataTable } from "../mocks/datatable";
import Collapsible from 'react-collapsible';
import Accordion from 'react-bootstrap/Accordion'
/*
Three tracks, exons, domains, secondary structres and disorder, exonTrack
*/
class NightingaleExonsVista extends Component {
  constructor(props){
    super(props);
    //console.log ('props',this.props)
    this.state = {
      sequence: this.props.sequence,
      domains: this.props.nightDom,
      ss: this.props.nightSs,
      exons: this.props.nighthTable,
  }}
  componentDidMount() {
     console.log("cdm1",this.props.transwhole)
     console.log("cdm2",this.props.transwhole.nightDom)
     //this["interpro-track"].innerHTML = this.props.transwhole.nightDom;
     let seqTrack="#sequence-track" + this.props.variable
     console.log(seqTrack)
     document.querySelector("#sequence-track"+this.props.variable).data = this.props.transwhole.sequence;
     document.querySelector("#domain-track"+this.props.variable).data = this.props.transwhole.nightDom;
     document.querySelector("#exon-track"+this.props.variable).data = this.props.transwhole.nightExons;
     document.querySelector("#secondaryStructure-track"+this.props.variable).data = this.props.transwhole.nightSs;
     document.querySelector("#disorder-track"+this.props.variable).data = this.props.transwhole.nightDis;     

    // document.querySelector("#saver2"+this.props.variable).backgroundColor = "#ddddee";
    }

    render () {
    loadWebComponent("protvista-manager", ProtvistaManager);
    loadWebComponent("protvista-track", ProtvistaTrack);
    loadWebComponent("protvista-navigation", ProtvistaNavigation);
    loadWebComponent("protvista-sequence", ProtvistaSequence);
    loadWebComponent("protvista-interpro-track", ProtvistaInterproTrack);
    loadWebComponent("protvista-saver", ProtvistaSaver);
    loadWebComponent("protvista-overlay", ProtvistaOverlay);
    loadWebComponent("protvista-zoom-tool", ProtvistaZoomTool);
    loadWebComponent("protvista-datatable", ProtvistaDatatable);
    console.log("render",this.props.nightDom)
    console.log("render2",this.props)
    console.log("STATE",this.state)
    console.log("TABLE", this.props.transwhole.nightExons)
    const lengthSeq=this.props.transwhole.sequence.length
    return (
      <div className="container" style = {{borderStyle:"double", padding: "0 0 20px 0",  borderWidth: "0 0 5px 0" /* top right bottom left */ }}>
          <Fragment>
            {/* <protvista-saver
              element-id="just-tracks" id={"saver2"+this.props.variable} file-name="tracks" file-format="jpeg">
              <button>Download Tracks</button>
            </protvista-saver> */}
            <protvista-manager
              displaystart= "1" displayend="100" length ={lengthSeq}>
              <div style={{padding:"10px", margin:"10px", paddingBottom:"20px"}}>
                <div style={{float:"left",
                  background: "#00639a",
                  color: "white",
                  background: "#00a6d5",
                  borderRadius: "4px",
                  margin: "2px"
                 }}>
                <span>{this.props.transwhole.tId}</span>
                </div>
              
              <protvista-zoom-tool
                length={lengthSeq}
                displaystart="1"
                displayend="100"
                style={{
                  float: "right",
                  "--button-background": "#00639a",
                  "--button-text-color": "#FFFFFF",
                  "--button-background-focus": "#00a6d5",
                  "--button-border-radius": "4px"}}>
                  <span slot="zoom-in">+</span>
                  <span slot="zoom-in-seq">Zoom to Sequence</span>
              </protvista-zoom-tool>
              </div>
              
              <div style={{padding:"20px"}}>
              <protvista-navigation length ={lengthSeq}  /> 
                  <div id={"just-tracks"}>
                      <protvista-sequence
                        length={lengthSeq} displaystart="1" displayend="100"
                        id={"sequence-track"+ this.props.variable}
                        highlight-event="onmouseover"
                        use-ctrl-to-zoom />
                      <protvista-track
                          id={"exon-track"+ this.props.variable}
                          length={lengthSeq}
                          // layout="non-overlapping"
                          highlight-event="onmouseover"
                          expanded
                          use-ctrl-to-zoom
                        />
                      <protvista-track
                          id={"secondaryStructure-track" + this.props.variable}
                          length={lengthSeq}
                          layout="non-overlapping"
                          use-ctrl-to-zoom
                        />
                        <protvista-track
                          id={"disorder-track" + this.props.variable}
                          length={lengthSeq}
                          layout="non-overlapping"
                          use-ctrl-to-zoom
                        />
                        <protvista-interpro-track
                          id={"domain-track" + this.props.variable}
                          length={lengthSeq}
                          highlight-event="onmouseover"
                          expanded
                          use-ctrl-to-zoom
                        />
                    </div>
                </div>
                <Accordion>
                  <Accordion.Item eventKey="0">
                    <Accordion.Header bsPrefix = "h5">
                      <div className="p-0 m-0 border">
                      Show/Hide Details
                      </div>
                    </Accordion.Header>
                  <Accordion.Body>
                    {/* <h5 style={{paddingBottom:"20px"}}>DataTable</h5> */}
                  <protvista-datatable displaystart="1" displayend="100">
                    <table>
                      <thead>
                        <tr>
                          <th data-filter="ft_key">Feature key </th>
                          <th>Description </th>
                          <th>Positions </th>
                        </tr>
                      </thead>
                      <tbody>
                        {/* {console.log("dataNew",dataTable)} */}
                        {this.props.transwhole.nighthTable?.map((row, i) => (
                          //i is row index, starting from 0
                          <Fragment key={i}>
                            {console.log("da*",row,i)}
                            <tr
                              data-id={`${row.start}-${row.end}`}
                              data-start={row.start}
                              data-end={row.end}
                            >
                              {/* inan loop and row content is decalred */}
                              <td data-filter="ft_key" data-filter-value={row.type}>
                                {row.type}
                              
                              </td>
                              <td>{row.description}</td>
                              <td>
                                {row.start} - {row.end}
                              </td>
                            </tr>
                            <tr data-group-for={`${row.start}-${row.end}`}>
                              <td>{row.tooltip}</td>
                            </tr>
                          </Fragment>
                        ))}
                      </tbody>
                    </table>
                    
                  </protvista-datatable>
                  </Accordion.Body>
                  </Accordion.Item>
                  </Accordion>
                  
                
        </protvista-manager>
      </Fragment>
      </div>
    )
  }
  }

export default NightingaleExonsVista;
