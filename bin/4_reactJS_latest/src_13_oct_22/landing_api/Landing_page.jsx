import React, { Component } from "react";
import Form from "../forms/Form";
import AppfetchOb from "./API_act_fetch";
import AppApi from "./API_handler";

class AppLand extends Component {
  constructor(props) {
    super(props);
    this.state = {
      post: {
        name: ""
      },
      jobs: [],
      showForm: true
    };
  }
  //have 3 state variables, showForm will render the form, post(name) and jobs [] unknown

  handleChange = e => {
    const { name, value } = e.target;
    console.log(name,value,'nameValue')
    this.setState(prevState => ({
      post: {[name]: value },
    }));
    
  };
  //When a value will be added to the user form, this change function will be called and then new value will be added to prevstate ('' in default firsttime), e.g. [name]: value is, store the value in it.
  handleSubmit = e => {
    e.preventDefault();
    this.setState(prevState => ({
      jobs: [prevState.post],
      showForm: false
    }));
  };

  render() {
    //const urlfetch = "http://172.16.3.146/nextrap/";
    //const urlfetch = "http://localhost:8001";
    //const urlfetch = "http://14.139.227.206/nextrap";
    const urlfetch = "http://localhost:8000";
    console.log("landingpage", urlfetch, this.state);
    if (this.state.showForm) {
      return (
        <div className="d-flex h-100">
          <div className="w-50 h-100 border ">
            {/* this is not working */}
            <div >
              <Form
                handleChange={this.handleChange}
                post={this.state.post}
                handleSubmit={this.handleSubmit}
              />
            </div>
          </div>
          <div className="w-50 h-100 border">
            <div className="position-relative" align="center">
              <h3>will add image carousla here with hyperlink to about nomenclaure and turtrials</h3>
            </div>
          </div>
        </div>
      );
    } else {
      if (isNaN(this.state.jobs[0].name)) {
        console.log(
          "insidestring funct",
          this.state.jobs[0],
          isNaN(this.state.jobs[0])
        );
        return (
          <div className="post-container">
            {this.state.jobs.map((job, index) => (
              <div key={index}>
                <AppApi
                  key={index}
                  urlprefix={urlfetch}
                  urlfetch={`${urlfetch}/name/${job.name}/`}
                />
              </div>
            ))}
          </div>
        );
      } else {
        console.log("insidenumeric funct");
        return (
          <div className="post-container">
            {this.state.jobs.map((job, index) => (
              <div key={index}>
                <AppfetchOb
                  key={index}
                  urlfetch={`${urlfetch}/ncbid/${job.name}/`}
                />
                ;
              </div>
            ))}
          </div>
        );
      }
    }
  }
}

// keep enetring text information to the form, once it is added, and you are satisfied, click submit, the information added to the name varibel then will be pushed to the jobs stack and then the showForm state will be locked.
// once we are out of reneding the stateform page, follwing infomratino iwll be iomporant 
// if added ifdntifier will be alphabetical then , isNaN(this.state.jobs[0].name) will be true and will go to first rebnederd block, and iwll try to fetch url with name job.name <component will be 'AppApi'
// otherwse iof numerical then the else block , delaing here can be tricky as if url isnt fetched and other details, this componenet may go awry. coponent will be 'AppfetchOb'

export default AppLand;
