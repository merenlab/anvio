# Purpose

The purpose of this directory is to keep track of people who are a
part of the anvi'o project. If you see someone who is missing,
please add them into the relevant file.

# How to add a new individual?

If you are adding a new anvi'o developer or contributor into either
of these YAML files, please send a PR with the relevant changes,
which include,

* **Edit [DEVELOPERS.yaml](https://github.com/merenlab/anvio/blob/master/anvio/data/misc/PEOPLE/DEVELOPERS.yaml) _or_ [CONTRIBUTORS.yaml](https://github.com/merenlab/anvio/blob/master/anvio/data/misc/PEOPLE/DEVELOPERS.yaml) files to add
a new entry** for the new person (please benefit from previous examples).
* **Add a new photo** under the `AVATARS` directory. It should be
a 900px x 900px head-shot (see previous examples).

Most fields are optional, but please note that every person mentioned
in either YAML files must have a GitHub username.

## A template to make things easeir

Feel free to use this **template** for new entries:

```
- github: (github username)
  name: (full name)
  twitter: (twitter username, if you have one and/or wish to list one)
  web: (https://your-web-page)
  avatar: (your-avatar.png that is in the AVATARS directory)
  email: (email address)
  linkedin: (linkedin username -- if there is one)
  orcid: (ORCiD, not the url, just the numbers -- if there is one)
  bio: "A one-sentence bio -- see examples at https://anvio.org/people/"
  affiliations:
    - title: (your title: Graduate Student / Post-doctoral scientist / Assistant Professor / Group Leader / etc)
      inst: (name of the institution you are affiliated with)
      inst_link: (the link to the institution or group page)
      current: (if this is a 'current' affiliation, put the word true here, if not, remove the line completely)
```

The purpose of affiliations is to keep track of every anvi'o person
starting from their first contribution to the platform. You can add
as many affiliations as you like. Please order them in such a way
that they are ordered from new to old.

# General things to consider developing this resource further

## Developers vs Contributors

`DEVELOPERS.yaml` contains individuals who has made literal contributions
to the anvi'o codebase. The information in this file can be used to
tag people in `__authors__` directives in anvi'o programs under the
[bin/](https://github.com/merenlab/anvio/tree/master/bin) and
[sandbox/](https://github.com/merenlab/anvio/tree/master/sandbox)
directories, which then can be used via the [anvio/authors.py](https://github.com/merenlab/anvio/blob/master/anvio/authors.py)
module that serves other programs such as [anvi-script-gen-help-pages](https://github.com/merenlab/anvio/blob/master/sandbox/anvi-script-gen-help-pages)
that generates anvi'o help pages at https://anvio.org/help

`CONTRIBUTORS.yaml` contains individuals has made indirect contributions
to the anvi'o community, including writing blog posts or tutorials, or
pushing the boundaries of the platform with their intellectual
contributions and guidance.

## Intellectual contributions

We have manually curated 'intellectual contributions' sections for
individuals listed in both [DEVELOPERS.yaml](https://github.com/merenlab/anvio/blob/master/anvio/data/misc/PEOPLE/DEVELOPERS.yaml) and [CONTRIBUTORS.yaml](https://github.com/merenlab/anvio/blob/master/anvio/data/misc/PEOPLE/DEVELOPERS.yaml).
At any given time, they will be inevitably incomplete. So please help
expand these entries to keep track of key contributions even for
those who are no longer around the project. Here are a few right-hand
rules to add more:

* Keep each contributions statement as concise and accurate as possible.
HTML notation is supported, but please limit the use of HTML to `a`
tags to link relevant documentation, PRs, help pages, or articles from
within the contrib statement.

* Each contributions statement must list a single contribution.

* If someone implemented an entirely new concept, program, workflow,
document, or framework in anvi'o, use the term 'Spearheaded'. I.e.,
"Spearheaded the development of anvi'o snakemake workflows". This
indicates that others may have contributed to it, or may contribute
in the future.

* If someone made notable contributions anything listed above, use
the phares "Made significant contributions". I.e., "Made significant
contributions to anvi'o snakemake workflows". What is notable and
what is not notable is not clear. Common sense is your best guidance,
but when you are unsure, discuss a given case with other developers.
A contribution may be notable only by you, so if you see something,
say something!
