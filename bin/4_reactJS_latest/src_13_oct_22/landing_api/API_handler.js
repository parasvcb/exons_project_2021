import React, { Component } from "react";
import AppfetchOb from "./API_act_fetch";
class AppApi extends Component {
  constructor(props) {
    super(props);
    this.state = {
      gene: false,
      showGeneListTable: false,
      url: "http://localhost:8000",
      genestate: ""
    };
    this.handleClick = this.handleClick.bind(this);
  }

  handleClick(genecode) {
    this.setState({ genestate: genecode, showGeneListTable: !this.state.showGeneListTable });
    console.log("clicked", this.state.showGeneListTable, genecode);
  }

  //remove constructor process at the top
  async componentDidMount() {
    try {
      const res = await fetch(`${this.props.urlfetch}`);
      console.log("fetchionglist url was");
      const todos = await res.json();
      console.log("fetchionglist", res, todos);
      this.setState({
        gene: todos,
        showGeneListTable: true
      });
    } catch (e) {
      console.log("e", e);
    }
  }
  /*
This module will only fetch the data from the api 
and then shift the control as well as data to 
next set of programs in spring folder, 
  componentDidMount()
    recieved props url will be rendered into JSON object, and added to state variable 'gene', showGeneListTable will be ON, 
    above state combination will reder the genetable, 
    
    clicking any entry from geneTable will switch on genestate variable of sttae and toggle the showGeneListTable,
    this combination will fire appfetch ob
    () details missingi n this module

*/
  //, <Apphd gene={this.state.gene} />
  render() {
    let testob = "";
    console.log("THIOSPROPS", this.props, this.state.showGeneListTable);
    if (this.state.gene && this.state.showGeneListTable) {
      if (this.state.gene.length === 0) {
        return (
          <h4>
            Sorry, there are no results for your query, please try again with
            another query
          </h4>
        );
      } else {
        return (
          <table className="table">
            <thead>
              <tr>
                <th>
                  <h5>Organism</h5>
                </th>
                <th>
                  <h5>Gene Id</h5>
                </th>
                <th>
                  <h5>Name</h5>
                </th>
                {/* <th>
                  <h5>Number of Transcripts</h5>
                </th> */}
                {/* <th>
                  <h5>Number of Exons</h5>
                </th> */}
              </tr>
            </thead>
            <tbody>
              {this.state.gene.map(gn => (
                <tr key={gn.id}>
                  <td>{gn.organism}</td>
                  <td>
                    <button
                      key={gn.id}
                      onClick={e => this.handleClick(gn.entrezid, e)}
                      className="m-2 btn btn-warning"
                    >
                      {gn.entrezid}
                    </button>
                  </td>
                  <td>{gn.name}</td>
                  {/* <td>{gn.trans.length}</td> */}
                  {/* <td>{gn.exonsGenes.length}</td> */}
                </tr>
              ))}
              ;
            </tbody>
          </table>
        );
      }
    } else {
      if (this.state.gene && this.state.genestate && !this.state.showGeneListTable) {
        console.log(
          "codes",
          this.state.genestate,
          this.state.gene[0].entrezid,
          this.state.gene
        );

        testob = this.state.gene.filter(mg => {
          return mg.entrezid === this.state.genestate;
        });
        return (
          <AppfetchOb
            urlfetch={`${this.props.urlprefix}/gene/${testob[0].id}/`}
          />
        );
      } else {
        return <p>Component loading....</p>;
      }
    }
  }
}
export default AppApi;
