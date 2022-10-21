import React, { Component } from "react";
import HolderFront from "../gene_render/HolderFront";

class AppfetchOb extends Component {
  state = {
    gene: false
  };
  //remove constructor process at the top
  async componentDidMount() {
    try {
      console.log("fresh", this.props)
      const res = await fetch(`${this.props.urlfetch}`);
      const todos = await res.json();

      this.setState({
        gene: todos
      });
    } catch (e) {
      console.log(e);
    }
  }

  /*
This module will only fetch the data from the api 
and then shift the control as well as data to 
next set of programs in spring folder
*/
  //, <Apphd gene={this.state.gene} />
  render() {
    console.log("adre to fethc", this.props);
    if (this.state.gene) {
      console.log("genestate", this.state.gene);
      if (
        "detail" in this.state.gene &&
        this.state.gene.detail === "Not found."
      ) {
        return (
          <h4>
            Sorry, there are no results for entered query, please try again with
            another query
          </h4>
        );
      } else {
        //console.log("geneob", this.state.gene);
        return <HolderFront gene={this.state.gene} />;
      }
    } else {
      return <p>Component fetching gene list... </p>
    }
  }
}

export default AppfetchOb;
