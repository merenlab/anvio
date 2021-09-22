

## Common errors

If you get an error that looks like this:
```
Config Error: The database at [PROFILE.db] does not seem to have a table named
              `detection_splits` :/ Here is a list of table names this database knows:
              [...]
```

That means your profile databases are not the correct version. The tables we are accessing in this program were introduced in profile database version 36. So the solution to this error is to migrate your databases to at least that version, using %{anvi-migrate}s. :) 
