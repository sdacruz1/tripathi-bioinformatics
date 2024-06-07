class Command {
    constructor(title, isEnabled, serverCall) {
        this.title = title;
        this.isEnabled = isEnabled;
        this.serverCall = serverCall; // Empty serverCall implies a section header
    }
}

class Parameter {
    constructor(title, isSub, type, options, serverTag, selected) {
        this.title = title;
        this.isSub = isSub;
        this.type = type;
        this.options = options;
        this.serverTag = serverTag;
        this.selected = selected;
    }
}

module.exports = {Command, Parameter};